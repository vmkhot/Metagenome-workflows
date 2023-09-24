library(DESeq2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(pheatmap)
library(gridExtra)
library(data.table)

setwd("./DESeq2/")

# read in column data
colData <- read.csv("./colData_samples_deseq2.csv",header = TRUE, sep = ',')
rownames(colData) <- colData[,1]

# factor column data by group and hours
factors <- c("condition","hours","diel","day","stage")
colData[factors] <- lapply(colData[factors], factor)

# read in count data - colData row names and countData column names have to match and be in exact order
countData <- read.csv("./countData_for_sample_plots.csv",header = TRUE, sep = ',')
rownames(countData) <- countData[,1]
countData[,1] <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create the DESeqDataSet object
ddsMat <- DESeqDataSetFromMatrix(countData, colData, design = ~ condition + day)

# run the differential expression analysis
dds <- DESeq(ddsMat, fitType='local')
resultsNames(dds)

sizeFactors(dds)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vsd <- vst(dds, blind = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA

pcaData <- plotPCA(vsd, intgroup=c("condition", "stage"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=stage, shape=condition)) +
  geom_point(size=5) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_manual(values=c("#3A9F1F","purple4","#226FD4"))

#"#2E8364","#287587","#2167A9",
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$day, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
Heatmap(sampleDistMatrix,
        clustering_distance_rows=sampleDists,
        clustering_distance_columns=sampleDists,
        row_split = 3, column_split = 3,
        col=colors)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

