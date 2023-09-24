library(ggplot2)

setwd("./DESeq2/cyano/")
df <- read.csv("./DE_percent_samples.csv",header = TRUE, sep = ',')


df$Time <- as.numeric(df$Time)

ggplot(df,aes(x = Time, y = isDE, group=1)) + geom_line(size=1.5, colour="royalblue") + geom_point(size=3,colour="darkblue") +
  labs(x = "Hours", y = "% of Genes Differentially Expressed")  +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank())+
  scale_x_continuous(breaks=df$Time, labels=df$Time)+
  annotate("rect", xmin=11,xmax=23,ymin=0,ymax=40,alpha=0.2, fill="darkblue")+
  annotate("rect", xmin=35,xmax=47,ymin=0,ymax=40,alpha=0.2, fill="darkblue")+
  annotate("rect", xmin=59,xmax=71,ymin=0,ymax=40,alpha=0.2, fill="darkblue")+
  annotate("rect", xmin=83,xmax=85,ymin=0,ymax=40,alpha=0.2, fill="darkblue")
  
  
  geom_rect(aes(df$diel))
            
  annotate("rect", xmin=12,xmax=18,ymin=0,ymax=39,alpha=0.2)
