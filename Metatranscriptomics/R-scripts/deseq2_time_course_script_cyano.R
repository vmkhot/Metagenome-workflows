library(DESeq2)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(gridExtra)
library(data.table)
library(tidyverse)
library("genefilter")

setwd("./DESeq2/cyano/")

# read in column data
colData <- read.csv("./colData_samples_deseq2_mod.csv",header = TRUE, sep = ',')
rownames(colData) <- colData[,1]

# factor column data by group and hours
factors <- c("condition","hours")
colData[factors] <- lapply(colData[factors], factor)

# read in count data - colData row names and countData column names have to match and be in exact order
countData <- read.csv("./countData_deseq2_mod.csv",header = TRUE, sep = ',')
rownames(countData) <- countData[,1]
countData[,1] <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DESEQ2
# create the DESeqDataSet object
ddsMat <- DESeqDataSetFromMatrix(countData, colData, design = ~ condition + hours + condition:hours)
# run the differential expression analysis
ddsTC <- DESeq(ddsMat, test="LRT", reduced = ~ condition + hours)

# parameters and query names
resultsNames(ddsTC)
sizeFactors(ddsTC)

# pull results from ddsTC (can query)
# althypothesis is that log2foldchange > lcfthreshold=0.5 (greatABS=two tailed)
# alpha is significance - has to match the padj later

set.seed(123)
res_0_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "condition_Treatment_vs_Control")
res_12_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours12")
res_24_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours24")
res_30_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours30")
res_36_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours36")
res_40_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours40")
res_44_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours44")
res_48_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours48")
res_52_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours52")
res_56_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours56")
res_60_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours60")
res_66_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours66")
res_68_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours68")
res_72_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours72")
res_78_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours78")
res_84_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours84")
res_96_ape <- lfcShrink(ddsTC, lfcThreshold=0.5, type = "apeglm", coef = "conditionTreatment.hours96")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot MA-plots and pull the data into a CSV
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
par(mfrow=c(4,4), mar=c(4,4,2,1))
ylim <- c(-3,3)

plotMA(res_0_ape,main = "res_0",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_12_ape,main = "res_12",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_24_ape,main = "res_24",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_30_ape,main = "res_30",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_36_ape,main = "res_36",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_40_ape,main = "res_40",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_44_ape,main = "res_44",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_48_ape,main = "res_48",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_52_ape,main = "res_52",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_56_ape,main = "res_56",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_60_ape,main = "res_60",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_66_ape,main = "res_66",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_68_ape,main = "res_68",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_72_ape,main = "res_72",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_78_ape,main = "res_78",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_84_ape,main = "res_84",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()
plotMA(res_96_ape,main = "res_96",alpha = 0.005, ylim=ylim,cex=0.7,colNonSig = "gray70",colSig = "royalblue"); drawLines()

res_0_df <- plotMA(res_0_ape,main = "res_0",returnData=TRUE)
write.csv(res_0_df, "res0_data_plotMA.csv")

res_12_df <- plotMA(res_12_ape,main = "res_12",returnData=TRUE)
res_24_df <- plotMA(res_24_ape,main = "res_24",returnData=TRUE)
res_30_df <- plotMA(res_30_ape,main = "res_30",returnData=TRUE)
res_36_df <- plotMA(res_36_ape,main = "res_36",returnData=TRUE)
res_40_df <- plotMA(res_40_ape,main = "res_40",returnData=TRUE)
res_44_df <- plotMA(res_44_ape,main = "res_44",returnData=TRUE)
res_48_df <- plotMA(res_48_ape,main = "res_48",returnData=TRUE)
res_52_df <- plotMA(res_52_ape,main = "res_52",returnData=TRUE)
res_56_df <- plotMA(res_56_ape,main = "res_56",returnData=TRUE)
res_60_df <- plotMA(res_60_ape,main = "res_60",returnData=TRUE)
res_66_df <- plotMA(res_66_ape,main = "res_66",returnData=TRUE)
res_68_df <- plotMA(res_68_ape,main = "res_68",returnData=TRUE)
res_72_df <- plotMA(res_72_ape,main = "res_72",returnData=TRUE)
res_78_df <- plotMA(res_78_ape,main = "res_78",returnData=TRUE)
res_84_df <- plotMA(res_84_ape,main = "res_84",returnData=TRUE)

dfList <- list(res_12_df,
               res_24_df,
               res_30_df,
               res_36_df,
               res_40_df,
               res_44_df,
               res_48_df,
               res_52_df,
               res_56_df,
               res_60_df,
               res_66_df,
               res_68_df,
               res_72_df,
               res_78_df,
               res_84_df)

dfList_new <- lapply(dfList, function(x) {
  df_out <- cbind(rownames(x), data.frame(x), row.names=NULL)
  df_out})

dfList_new %>%
  imap(function(x, y) x %>% rename_with(~paste(., y, sep = '_'), -"rownames(x)")) %>% 
  reduce(full_join, by = 'rownames(x)') -> final_df

write.csv(final_df, "merged_data_plotMA.csv",sep = ",")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# library("genefilter")
vsd <- assay(vst(ddsTC))
topVarGenes <- head(order(rowVars(vsd)), decreasing = TRUE, 200)
mat  <- vsd[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
Heatmap(mat, cluster_columns = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PULL CRISPR Genes significant in res_44_ape

# extract log2 fold changes using coef()
betas <- coef(ddsTC)
betas_df <- as.data.frame(betas)
colnames(betas)
#betas_f <- cbind(rownames(betas), data.frame(betas), row.names=NULL)

# filter betas dataframe by svalue < 0.005 in res_44 & res_72 samples
betas_filtered <- betas_df %>% 
  filter(res_44_ape$svalue < 0.005)
betas_filtered <- betas_filtered[,19:33]

#filter out just CRISPRs
include_list <- c("CRISPR", "crispr_repeat")
betas_CRISPR <- betas_filtered %>% 
  filter(rownames(betas_filtered) %like% ("CRISPR|crispr"))

Heatmap(betas_CRISPR,
        cluster_columns = FALSE)

betas_toxin <- betas_filtered %>% 
  filter(rownames(betas_filtered) %like% ("toxin"))
Heatmap(betas_toxin,
        cluster_columns = FALSE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PULL normalized counts for clustering samples
library(limma)
vsd <- vst(ddsTC)
mat <- assay(vsd)
mm <- model.matrix(~ condition + hours, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
plotPCA(vsd)
