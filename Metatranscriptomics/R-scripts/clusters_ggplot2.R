library(tidyverse)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)

setwd("./DESeq2/cyano/clustering_genes/all_genes/")

all_data <- read.csv("replotting_clusters.csv",sep=",", header=TRUE)

colnames(all_data) <- c("geneid","cluster","0","12","24","30","36","40","44","48","52","56","60","66","68","72","78","84")

all_data <- as_tibble(all_data)


subset_m <- all_data %>%
  filter(grepl("7", cluster, fixed = TRUE)) %>%
  replace(is.na(.),0) %>%
  pivot_longer(c("0","12","24","30","36","40","44","48","52","56","60","66","68","72","78","84"),
               names_to = "hour", values_to = "log2FC") 

subset_m_means <- subset_m %>%
  group_by(cluster,hour) %>%
  mutate(means = mean(log2FC))


subset_m_means$hour <- as.numeric(subset_m_means$hour)


p_sub <- ggplot(subset_m_means, aes(x=hour, y=log2FC, group=geneid))+
  #geom_ribbon(data=subset_m_means,aes(hour,ymin = means -2, ymax = means +2),fill="lightblue",alpha=0.25) +
  geom_hline(aes(yintercept = 0),linetype="dashed", lwd=0.75) +
  geom_line(size=0.75,color="grey40", alpha=0.8)+
  geom_smooth(data=subset_m_means,aes(hour,means),color="royalblue", se=FALSE,size=2)+
  geom_vline(aes(xintercept = 36),linetype="dotdash", lwd=0.75)+
  geom_vline(aes(xintercept = 66),linetype="dotdash", lwd=0.75)+
  scale_x_continuous(breaks=subset_m_means$hour, labels=subset_m_means$hour,guide = guide_axis(n.dodge = 2))+
  facet_wrap(~subset_m_means$cluster,
             labeller = label_wrap_gen(width = 20), nrow = 3)+
  ylim(-10,10)+
  labs(x="Hours", y="Log2 Fold Change")+
  theme(axis.text.y = element_text(colour = "black", size = 12,face="bold"), 
        axis.text.x = element_text(colour = "black", size = 12), 
        axis.title = element_text(face = "bold", size = 14, colour = "black"), 
        legend.position = "none", 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        strip.background = element_blank(),strip.text = element_text(colour = "black", size = 12, face="bold"))

p_sub











all_data_m <- all_data %>%
  replace(is.na(.),0) %>%
  pivot_longer(c("0","12","24","30","36","40","44","48","52","56","60","66","68","72","78","84"),
               names_to = "hour", values_to = "log2FC") %>%
  group_by(cluster,hour) %>%
  mutate(means = mean(log2FC))


all_data_m$hour <- as.numeric(all_data_m$hour)

p1 <- ggplot(all_data_m, aes(x=hour, y=log2FC, group=geneid))+
  geom_hline(aes(yintercept = 0),linetype="dashed", lwd=0.75) +
  geom_line(size=0.75,color="grey40", alpha=0.8)+
  geom_smooth(data=all_data_m,aes(hour,means),color="royalblue", se=FALSE,size=2)+
  geom_vline(aes(xintercept = 36),linetype="dotdash", lwd=0.75)+
  geom_vline(aes(xintercept = 66),linetype="dotdash", lwd=0.75)+
  scale_x_continuous(breaks=all_data_m$hour, labels=all_data_m$hour,guide = guide_axis(n.dodge = 2))+
  facet_wrap(~all_data_m$cluster,
             labeller = label_wrap_gen(width = 20), nrow = 3)+
  labs(x="Hours", y="Log2 Fold Change")+
  ylim(-10,10)+
  theme(axis.text.y = element_text(colour = "black", size = 12,face="bold"), 
        axis.text.x = element_text(colour = "black", size = 12), 
        axis.title = element_text(face = "bold", size = 14, colour = "black"), 
        legend.position = "none", 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        strip.background = element_blank(),strip.text = element_text(colour = "black", size = 12, face="bold"))

p1
library(svglite)
ggsave("clusters_ggplot2.svg", plot=p1, width = 500, height = 300,units = "px",rast)
