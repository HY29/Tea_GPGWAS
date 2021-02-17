############################################################
### GWAS flow 
### Draft script 200811 HY
### Update 210217 HY
############################################################
##load packages
library(rrBLUP)
library(qqman)
#set directory
setwd("/Users/HirotoYAMASHITA/BioInfo/Ranalysis/01_Genetics/TeaDemo/")

###Prepare dataframe----
##read phenotype
pheno <- read.csv("Pheno/PhenoAllMean.csv", row.names = 1, fileEncoding="CP932")
## read genotype
geno <- read.csv("Geno/list_gwas_maf0.05_miss0.7.imp.csv")
##data arrangement
chr <- geno$chr
pos <- geno$pos
id <- geno$id
gt.score <- geno[,5:154] #150 accessions
gt.score.t <- t(gt.score)
dim(gt.score.t)
#extract geno of accessions with pheno
gt.score.t <- gt.score.t[rownames(pheno),]
dim(gt.score.t)

##estimation genetic background by PCA 
pca <- prcomp(gt.score.t)
summary(pca)
##extraction cum
pca_sum <- summary(pca)
pca_cum <- as.data.frame(t(pca_sum$importance))
dir.create(path="GenoPCA")
write.csv(pca_cum, "GenoPCA/geno_pca_cum.csv")
pdf("GenoPCA/plot_pca_cum.pdf", width=3, height=3)
plot(pca)
dev.off()

##prepare genotype data for GWAS function
g <- data.frame(id, chr, pos, t(gt.score.t))
row.names(g) <- 1:nrow(g)
colnames(g) <- c("id", "chr", "pos", rownames(gt.score.t))
##amat
amat <- A.mat(gt.score.t, shrink=TRUE)
colnames(amat) <- row.names(amat) <- row.names(gt.score.t)
amat[1:6, 1:6]

###GWAS_ManhattanPlot_QQplot----
#make directory
dir.create(path="GWAS")
dir.create(path="GWAS/MLM_QK_pc6")
dir.create(path="GWAS/MLM_QK_pc6/Fig")
dir.create(path="GWAS/MLM_QK_pc6/Table")
##select a trait (actually, you can select multiple traits here)
for (i in 1:length(pheno)){
  print(i)
  y = pheno[,i]
  Trait = colnames(pheno[i])
  ##prepare phenotypic values for GWAS
  p <- data.frame(rownames(gt.score.t), y)
  colnames(p) <- c("gid", "y")
  ## GWAS with Q (population structure) + K (kinship relatedness)
  ## MAF is set as 0.05
  gwa <- GWAS(p, g, K=amat, n.PC=6, min.MAF = 0.05, plot = F) #MLM
  #gwa <- GWAS(p, g, min.MAF = 0.05, plot = F) #GLM
  
  ##replace -log P=0 with NA
  gwa$y[gwa$y==0] <-NA
  
  ##Data arrange for draw manhattan plot 
  mht <- data.frame(SNP=gwa$id, CHR=gwa$chr, BP=gwa$pos, P=10^(-gwa$y))
  mht <- na.omit(mht)
  main=paste("",Trait,"",sep="")
  
  #draw manhattan plot
  png(paste("GWAS/MLM_QK_pc6/Fig/pheno_",i,"_",Trait,"_manhattan.png", sep=""),width=350, height=300)
  manhattan(mht, 
          cex.lab=1.2,cex.axis=1.2, cex.main=1.5,cex=1.3,
          suggestiveline=FALSE,genomewideline=FALSE, main=main, col=c("dodgerblue4", "gray66"))
  dev.off()
  
  #output gwas results
  write.csv(mht,
          paste("GWAS/MLM_QK_pc6/Table/pheno_",i,"_",Trait,"_gwas_result.csv", sep=""))
  
  ##draw Q-Q plot
  pdf(file =paste("GWAS/MLM_QK_pc6/Fig/pheno_",i,"_",Trait,"_QQplot.pdf",sep=""), width=4, height=4)
  qq(mht$P, main=main,col="black", cex.lab=1,cex.axis=1,cex.main=1,cex=1)
  dev.off()
  
  }

###MultipleTest----
#calculate false discovery rate (BH method)
fdr <- p.adjust(mht$P, method="BH")
#select markers under FDR
fdr.thresh = 0.05

sig <- mht[fdr < fdr.thresh, ]
dim(sig)

# calculate false discovery rate (BH method)
#draw manhattan plot with coloring significant SNPs
manhattan(mht, highlight=sig$SNP)

#local manhattan
chr1 <- subset(mht, CHR==1)
manhattan(chr1, suggestiveline=F, col="turquoise4")

###LocalManhattanPlot----
#close up the chromosomes (for example chrom 1, 2)
manhattan(subset(mht, CHR %in% c(1,2), highlight = sig$SNP))

# calculate false discovery rate (Bonferroni method) 
p.bonferroni <- p.adjust(mht$P, method="bonferroni")
#select markers under the thresould
p.thresh=0.05

sig <- mht[p.bonferroni < p.thresh, ]
dim(sig)

#draw manhattan plot with coloring significant SNPs
manhattan(mht, highlight=sig$SNP)






