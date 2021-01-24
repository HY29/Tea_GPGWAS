#########################################################
### Genomic prediction by rrBLUP(RR)
### Draft script 200811 HY
#########################################################
library(rrBLUP)

###Prepare dataframe----
#read phenotype
pheno <- read.csv("pheno/phenotype_all_mean.csv", row.names = 1, fileEncoding="CP932")
#read genotype
geno <- read.csv("geno/list_gwas_maf0.05_miss0.7.imp.csv")
#data arrangement
gt.score <- geno[,5:154] #150 accessions
x <- t(gt.score)


###GenomicPrediction----
# set n-fold cross varidation
n.fold = 10
# set number of repeat
rep = 10
#make directory for GBLUP(GAUSS)
dir.create(path="GP/GBLUP(RR)")
model = print(paste("GBLUP(RR)")) ##Caution!

#Repeat this analysis for n.folds
for (z in 1:length(pheno)){
  print(z)
  y = pheno[,z]
  dir.create(path=paste("GP/",model,"/pheno_",z,"", sep=""))
for (r in 1:rep){
  print(paste(z,r))		
  id <- sample(1:length(y) %% n.fold)
  id[id == 0] <- n.fold
  
  y.pred <- rep(NA, times=length(y))
  for(i in 1:n.fold) {
  print(paste(z,r,i))
  y.train <- y[id != i]
  x.train <- x[id != i,]
  x.test <- x[id == i,]
  
  res <- kinship.BLUP(y.train, x.train, x.test, K.method = "RR")
  y.pred[id == i] <- res$g.pred + rep(res$beta, length(res$g.pred))
  }
  
  write.csv(y.pred, paste("GP/",model,"/pheno_",z,"/pheno_",z,"_GP_10-fold_gwSNP_", r,".csv", sep=""))
  
}}
