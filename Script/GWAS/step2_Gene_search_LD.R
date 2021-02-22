#########################################################
### Search of Candidate genes
### Draft script 200811 HY
### Update 210217 HY
#########################################################
#load package
library(tidyverse)
#set directory
setwd("/Users/HirotoYAMASHITA/BioInfo/Ranalysis/01_Genetics/TeaDemo/")

###10kb gene_GLM&MLM for single GWAS res----
#GeneInf
data<-read.csv("GeneInf/CSS_ChrLev_GeneInf_withoutContig.csv",header=T)

#TopSNPs by GWAS
pheno = paste("Ala") #Select phenotype
phenoNo = paste("1") #Select phenotype
#make directory
dir.create(path="GWAS/MLM_QK_pc6/CandidateGenes")

GWASresult <- read.csv(paste("GWAS/MLM_QK_pc6/Table/pheno_",phenoNo,"_",pheno,"_gwas_result.csv", sep=""), row.names = 1)
GWASresult <- GWASresult[order(GWASresult$P),] #sort pvalue
topsnp <- GWASresult[1:160,] ##Caution!! #Please select number of SNPs by GP curve
matassoSNPs <- data.frame(topsnp$CHR,topsnp$BP,topsnp$P)
colnames(matassoSNPs) <- c("Chr","position")

#Set of LD window and gene promoter region
region<-10000 #10kb #Please select the number based on LD decay
prom<-2000
end<-0

Fposi<-data[data$ori=="+",]
Fposi$up<-Fposi$up-prom
Fposi$down<-Fposi$down+end

Rposi<-data[data$ori=="-",]
Rposi$down<-Rposi$down+prom
Rposi$up<-Rposi$up-end

genemat<-rbind(Fposi,Rposi)
genemat<-genemat[order(genemat$Chr,genemat$up),]

#Emptymatrix for results
resmat<-data.frame()

#Extraction of associated genes by GWAS top-ranked SNPs
for (n in 1:nrow(matassoSNPs)){
  if(n%%1000==0){
    print(n)
  }
  SNPres<-data.frame()
  
  genematchr<-genemat[genemat$Chr==matassoSNPs$Chr[n],]
  
  res<-genematchr[genematchr$down>=(matassoSNPs[n,2]-region)&genematchr$up<=(matassoSNPs[n,2]+region),]
  resSNP<-rep(paste(matassoSNPs[n,1],"_",matassoSNPs[n,2],sep=""),nrow(res))
  resPval<-rep(matassoSNPs[n,3],nrow(res))
  dist<-rep(NA,nrow(res))
  
  
  selectw<-matassoSNPs[n,2]>=res$up&matassoSNPs[n,2]<=res$down
  selectuf<-matassoSNPs[n,2]<res$up&res$ori=="+"
  selectuR<-matassoSNPs[n,2]<res$up&res$ori=="-"
  selectdf<-matassoSNPs[n,2]>res$down&res$ori=="+"
  selectdR<-matassoSNPs[n,2]>res$down&res$ori=="-"
  
  dist[selectw]<-"within"
  dist[selectuf]<-(matassoSNPs[n,2]-res$up)[selectuf]
  dist[selectuR]<-(res$up-matassoSNPs[n,2])[selectuR]
  dist[selectdf]<-(matassoSNPs[n,2]-res$down)[selectdf]
  dist[selectdR]<-(res$down-matassoSNPs[n,2])[selectdR]
  
  SNPres<-cbind(resSNP,resPval,res,dist)
  
  resmat<-rbind(resmat,SNPres)
}

write.csv(resmat,
          paste("GWAS/MLM_QK_pc6/CandidateGenes/GWAS_Top",nrow(matassoSNPs),"SNPs_pheno_",phenoNo,"_",pheno,"_",region/1000,"kb_candidateGenes.csv",sep=""),
          quote=F,row.names=F)

###Combined to Genes annotation----
GenesAno<-read.csv("GeneInf/CSS_ChrLev_20200506_Function_nr.csv",header=T)
resmatAno <- dplyr::left_join(resmat, GenesAno, by="GeneID", na_matches="never")
write.csv(resmatAno,
          paste("GWAS/MLM_QK_pc6/CandidateGenes/GWAS_Top",nrow(matassoSNPs),"SNPs_pheno_",phenoNo,"_",pheno,"_",region/1000,"kb_candidateGenes_Anotation.csv",sep=""),
          quote=F,row.names=F)


