#########################################################
### Calculation of R and RMSE value in GPwithGWAS
### Draft script 200811 HY
### Update 210217 HY
#########################################################
#set directory
setwd("/Users/HirotoYAMASHITA/BioInfo/Ranalysis/01_Genetics/TeaDemo/")
#read phenotype
pheno <- read.csv("Pheno/PhenoAllMean.csv", row.names = 1, fileEncoding="CP932")

###Output of R & RMSE value in each round----
#set model
model = print(paste("500SNP_GBLUP(RR)")) ##Caution!

#all file read in directory
for (n in 1:length(nsnp)) {
  print(paste(n))
  for (z in 1:length(pheno)){
    print(paste(n,z))
    y = pheno[,z]
    files <- list.files(paste("GPwithGWAS/",model,"/pheno_",z,"/Top_",nsnp[n],"/", sep=""))
    
    mat<-matrix(NA,nrow=length(y),ncol=length(files)+1)
    for (r in 1:length(files)){
      print(paste(n,z,r))
      one<-read.csv(paste("GPwithGWAS/",model,"/pheno_",z,"/Top_",nsnp[n],"/GP_10-fold_topSNP_",r,".csv", sep=""),row.names=1)
      mat[,r]<-one[,1]
      names<-c(names,paste(r,sep=""))
      }
    #input actual pheno data to matrix of predicted data
    mat[,length(files)+1] <- y
    #R & RMSEvalue
    R.value<-c()
    RMSE.value<-c()
    for(m in 1:length(files)){
      print(paste(n,z,r,m))
      #Rvalue
      R<-cor(y,mat[,m])
      R.value<-c(R.value,R)
      print(R.value)
      #RMSEvalue
      RMSE<-sqrt((sum((y-mat[,m])^2))/length(y))
      RMSE.value<-c(RMSE.value,RMSE)
      print(RMSE.value)
      }
    write.csv(R.value, paste("GPwithGWAS/",model,"/pheno_",z,"/Top_",nsnp[n],"/GP_10-fold_topSNP_Rvalue.csv", sep=""))
    write.csv(RMSE.value, paste("GPwithGWAS/",model,"/pheno_",z,"/Top_",nsnp[n],"/GP_10-fold_topSNP_RMSEvalue.csv", sep=""))
    }}

###Calculation of mean of R & RMSE value----
#Mean & SD of R.value 
for (z in 1:length(pheno)) {
  print(paste(n))
  matR<-matrix(NA,nrow=length(R.value),ncol=length(nsnp))
  matR_MeanSD<-matrix(NA,ncol=2,nrow=length(nsnp))
  colnames(matR_MeanSD)<-c("R2_Mean","R2_SD")
  rownames(matR_MeanSD)<-nsnp
  for (n in 1:length(nsnp)){
    print(paste(z,n))
    one<-read.csv(paste("GPwithGWAS/",model,"/pheno_",z,"/Top_",nsnp[n],"/GP_10-fold_topSNP_Rvalue.csv", sep=""),row.names=1)
    matR[,n]<-one[,1]
    matR_MeanSD[n,1]<-mean(matR[,n])
    matR_MeanSD[n,2] <-sd(matR[,n])
    write.csv(matR_MeanSD, paste("GPwithGWAS/",model,"/pheno_",z,"/AllTopSNP_GP_10-fold_topSNP_MeanSD_Rvalue.csv", sep=""))
    }}

#Mean & SD of RMSE
for (z in 1:length(pheno)) {
  print(paste(n))
  matRMSE<-matrix(NA,nrow=length(R.value),ncol=length(nsnp))
  matRMSE_MeanSD<-matrix(NA,ncol=2,nrow=length(nsnp))
  colnames(matRMSE_MeanSD)<-c("RMSE_Mean","RMSE_SD")
  rownames(matRMSE_MeanSD)<-nsnp
  for (n in 1:length(nsnp)){
    print(paste(z,n))
    one<-read.csv(paste("GPwithGWAS/",model,"/pheno_",z,"/Top_",nsnp[n],"/GP_10-fold_topSNP_RMSEvalue.csv", sep=""),row.names=1)
    matRMSE[,n]<-one[,1]
    matRMSE_MeanSD[n,1]<-mean(matRMSE[,n])
    matRMSE_MeanSD[n,2] <-sd(matRMSE[,n])
    write.csv(matRMSE_MeanSD, paste("GPwithGWAS/",model,"/pheno_",z,"/AllTopSNP_GP_10-fold_topSNP_MeanSD_RMSEvalue.csv", sep=""))
    }}
