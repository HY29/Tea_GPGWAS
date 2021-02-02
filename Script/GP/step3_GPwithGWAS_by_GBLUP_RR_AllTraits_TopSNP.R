#########################################################
### GP modelling using top-ranked SNPs by GWAS
### Draft script 200811 HY
#########################################################
library(rrBLUP)


###Prepare dataframe----
#read phenotype
pheno <- read.csv("pheno/phenotype_all_mean.csv", row.names = 1, fileEncoding="CP932")
#read genotype
geno <- read.csv("geno/list_gwas_maf0.05_miss0.7.imp.csv")
#make directory for model
#dir.create(path="GPwithGWAS")
dir.create(path="GPwithGWAS/500SNP_GBLUP(RR)")
model = print(paste("500SNP_GBLUP(RR)")) ##Caution!

###GP using top-ranked SNPs by GWAS----
#Extraction top-ranked SNPs
for (z in 1:length(pheno)){
  print(z)
  y = pheno[,z]
  Trait = colnames(pheno[z])
  dir.create(path=paste("GPwithGWAS/",model,"/pheno_",z,"", sep=""))
  
  # read GWAS topsnp marker
  GWASresult <- read.csv(paste("GWAS/MLM_QK_pc6/pheno_",z,"_",Trait,"_gwas_result.csv", sep=""), row.names = 1)
  GWASresult <- GWASresult[order(GWASresult$P),] #sort pvalue
  topsnp <- GWASresult[1:500,] #extract of top30 results ### Caution !!! ###

  # extraction topsnp marke
  topsnp_geno <- geno[rownames(topsnp),] 
  dim(topsnp_geno)
  #data arrangement
  gt.score <- topsnp_geno[,5:154] #150 accessions
  xx <- t(gt.score)
  #set number of topsnp
  nsnp<-c(20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500)#Caution!!!
  
  #GP modelling
  for (n in 1:length(nsnp)) {
    print(paste(z,n))
    dir.create(path=paste("GPwithGWAS/",model,"/pheno_",z,"/Top_",nsnp[n],"", sep=""))
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
    
    write.csv(y.pred, paste("GPwithGWAS/",model,"/pheno_",z,"/Top_",nsnp[n],"/GP_10-fold_topSNP_",r,".csv", sep=""))
    
  }}}

