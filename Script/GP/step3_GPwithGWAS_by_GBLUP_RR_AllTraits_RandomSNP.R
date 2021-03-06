#########################################################
### GP modelling using top-ranked SNPs by GWAS
### Draft script 200811 HY
### Update 210217 HY
#########################################################
#load package
library(rrBLUP)
#set directory
setwd("/Users/HirotoYAMASHITA/BioInfo/Ranalysis/01_Genetics/TeaDemo/")

###Prepare dataframe----
#read phenotype
pheno <- read.csv("Pheno/PhenoAllMean.csv", row.names = 1, fileEncoding="CP932")
#read genotype
geno <- read.csv("Geno/list_gwas_maf0.05_miss0.7.imp.csv")

#select Model
model = print(paste("500SNP_GBLUP(RR)")) ##Caution!
#make directory for model
dir.create(path="GPwithGWAS")
dir.create(path=paste("GPwithGWAS/",model,"_Random", sep=""))

###GP using top-ranked SNPs by GWAS----
#Extraction top-ranked SNPs
set.seed(100)
for (z in 1:length(pheno)){
  print(z)
  y = pheno[,z]
  Trait = colnames(pheno[z])
  dir.create(path=paste("GPwithGWAS/",model,"_Random/pheno_",z,"", sep=""))
  
  # read GWAS topsnp marker
  set.seed(10)
  GWASresult <- read.csv(paste("GWAS/MLM_QK_pc6/Table/pheno_",z,"_",Trait,"_gwas_result.csv", sep=""), row.names = 1)
  GWASresult <- GWASresult[sample(rownames(GWASresult)),] #sort pvalue
  topsnp <- GWASresult[1:500,] #extract of top30 results ##Caution!!

  # extraction topsnp marke
  topsnp_geno <- geno[rownames(topsnp),] 
  dim(topsnp_geno)
  #data arrangement
  gt.score <- topsnp_geno[,5:154] #150 accessions
  xx <- t(gt.score)
  #set number of topsnp
  nsnp<-seq(20,500,by=20) #Caution!!!
  
  #GP modelling
  for (n in 1:length(nsnp)) {
    print(paste(z,n))
    dir.create(path=paste("GPwithGWAS/",model,"_Random/pheno_",z,"/Top_",nsnp[n],"", sep=""))
    x <-xx[, 1:nsnp[n]]

    # set n-fold cross varidation
    n.fold = 10
    # set number of repeat
    rep = 10
   
      for (r in 1:rep){
        print(paste(z,n,r))		
        id <- sample(1:length(y) %% n.fold)
        id[id == 0] <- n.fold
    
        y.pred <- rep(NA, times=length(y))
        for(i in 1:n.fold) {
          print(paste(z,n,r,i))
            y.train <- y[id != i]
            x.train <- x[id != i,]
            x.test <- x[id == i,]
      
            res <- kinship.BLUP(y.train, x.train, x.test, K.method = "RR")
            y.pred[id == i] <- res$g.pred + rep(res$beta, length(res$g.pred))
    }
    
    write.csv(y.pred, paste("GPwithGWAS/",model,"_Random/pheno_",z,"/Top_",nsnp[n],"/GP_10-fold_topSNP_",r,".csv", sep=""))
    
  }}}

