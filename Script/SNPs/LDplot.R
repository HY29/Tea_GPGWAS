############################################################
### LD decay estimation 
### Draft script 200811 HY
############################################################
#load packages
library(ggplot2)
library(cowplot)

#Data arrangement
d <- read.table("list_gwas_maf0.05_miss0.7_ld_window_50kb.geno.ld", header=TRUE)
distance <- as.data.frame(d$POS2-d$POS1)
LD <- data.frame(distance, d$R.2)
colnames(LD) <- c("distance","R2")
plot(LD)
write.csv(LD, "LD_miss0.7_w50kb.csv")

#Mean R2
mean(LD$R2)

#Visualization by ggplot2
pdf("LDplot_miss0.7.pdf",width=3, height=3)

ggplot(LD, aes(y=R2, x=distance))+
  #guides(color=FALSE)+
  geom_point(size=1)+
  geom_smooth(method="loess")+
  theme_cowplot(font_size = 12, line_size = 1.0)+
  theme(axis.text.x = element_text(angle=0))+
  ylab("Linkage disequilibrium (r2)")+
  xlab("Distance (bp)")

dev.off()
