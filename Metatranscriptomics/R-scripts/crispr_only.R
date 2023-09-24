library(tidyverse)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)

setwd("../../")

all_data <- read.csv("crispr_only_genes.csv",sep=",", header=TRUE)

colnames(all_data) <- c("geneid","cluster","name","regulation","dynamic","0","12","24","30","36","40","44","48","52","56","60","66","68","72","78","84")

all_data <- as_tibble(all_data)

all_data_m <- all_data %>%
  filter(!grepl("-1", regulation, fixed = TRUE)) %>%
  pivot_longer(c("0","12","24","30","36","40","44","48","52","56","60","66","68","72","78","84"),
               names_to = "hour", values_to = "log2FC")

all_data_m$hour <- as.numeric(all_data_m$hour)

p1 <- ggplot(all_data_m, aes(x=hour, y=log2FC, group=geneid), color=name)+
  geom_hline(aes(yintercept = 0),linetype="dashed", lwd=0.75) +
  geom_point(aes(colour=name),size=1.3,alpha=0.8)+
  geom_line(aes(colour=name),size=0.75,alpha=0.8)+
  geom_vline(aes(xintercept = 36),linetype="dotdash", lwd=0.75)+
  geom_vline(aes(xintercept = 66),linetype="dotdash", lwd=0.75)+
  scale_colour_brewer(breaks = all_data_m$name, palette = "Paired")+
  scale_x_continuous(breaks=all_data_m$hour, labels=all_data_m$hour,guide = guide_axis(n.dodge = 2))+
  facet_wrap(all_data_m$dynamic~all_data_m$regulation,
             labeller = label_wrap_gen(width = 20), nrow = 3, scales = "free")+
  labs(x="Hours", y="Log2 Fold Change",colour="CRISPR-Associated Gene")+
  ylim(-10,10)+
  theme(axis.text.y = element_text(colour = "black", size = 12,face="bold"), 
        axis.text.x = element_text(colour = "black", size = 12), 
        axis.title = element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        strip.background = element_blank(),strip.text = element_text(colour = "black", size = 12, face="bold"),
        legend.key = element_blank())

p1
