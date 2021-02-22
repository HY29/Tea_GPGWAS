#########################################################
### Calculation of R and RMSE value in genomic prediction
### Draft script 200811 HY
### Update 210217 HY
#########################################################
#set directory
setwd("/Users/HirotoYAMASHITA/BioInfo/Ranalysis/01_Genetics/TeaDemo/")
#read phenotype
pheno <- read.csv("Pheno/PhenoAllMean.csv", row.names = 1, fileEncoding="CP932")

###Output of R & RMSE value in each round----
#set model
model = print(paste("GBLUP(RR)"))

#all file read in directory
for (z in 1:length(pheno)){
    print(z)
    y = pheno[,z]
    files <- list.files(paste("GP/",model,"/pheno_",z,"/", sep=""))
    
    mat<-matrix(NA,nrow=length(y),ncol=length(files)+1)

for (r in 1:length(files)){
  print(paste(z,r))
  one<-read.csv(paste("GP/",model,"/pheno_",z,"/pheno_",z,"_GP_10-fold_gwSNP_", r,".csv", sep=""),row.names=1)
  mat[,r]<-one[,1]
  names<-c(names,paste(r,sep=""))
}

#input actual pheno data to matrix of predicted data
mat[,length(files)+1] <- y

#R & RMSEvalue
R.value<-c()
RMSE.value<-c()

for(m in 1:length(files)){
  print(paste(z,r,m))
  R<-cor(y,mat[,m])
  R.value<-c(R.value,R)
  print(R.value)

  RMSE<-sqrt((sum((y-mat[,m])^2))/length(y))
  RMSE.value<-c(RMSE.value,RMSE)
  print(RMSE.value)
}
  write.csv(R.value, paste("GP/",model,"/pheno_",z,"/pheno_",z,"_GP_10-fold_gwSNP_Rvalue.csv", sep=""))
  write.csv(RMSE.value, paste("GP/",model,"/pheno_",z,"/pheno_",z,"_GP_10-fold_gwSNP_RMSEvalue.csv", sep=""))
  
}

###Calculation of mean of R & RMSE value----
# Mean & SD of R.value
matR<-matrix(NA,nrow=length(R.value),ncol=length(pheno))
matR_MeanSD<-matrix(NA,ncol=2,nrow=length(pheno))
colnames(matR_MeanSD)<-c("R2_Mean","R2_SD")
rownames(matR_MeanSD)<-colnames(pheno)
for (i in 1:length(pheno)){
  print(paste(i))
  one<-read.csv(paste("GP/",model,"/pheno_",i,"/pheno_",i,"_GP_10-fold_gwSNP_Rvalue.csv", sep=""),row.names=1)
  matR[,i]<-one[,1]
  matR_MeanSD[i,1]<-mean(matR[,i])
  matR_MeanSD[i,2]<-sd(matR[,i])
  
}
  write.csv(matR_MeanSD, paste("GP/",model,"/Allpheno_GP_10-fold_gwSNP_MeanSD_Rvalue.csv", sep=""))

# Mean & SD of RMSE
matRMSE<-matrix(NA,nrow=length(R.value),ncol=length(pheno))
matRMSE_MeanSD<-matrix(NA,ncol=2,nrow=length(pheno))
colnames(matRMSE_MeanSD)<-c("RMSE_Mean","RMSE_SD")
rownames(matRMSE_MeanSD)<-colnames(pheno)
for (i in 1:length(pheno)){
  print(paste(i))
  one<-read.csv(paste("GP/",model,"/pheno_",i,"/pheno_",i,"_GP_10-fold_gwSNP_RMSEvalue.csv", sep=""),row.names=1)
  matRMSE[,i]<-one[,1]
  matRMSE_MeanSD[i,1]<-mean(matRMSE[,i])
  matRMSE_MeanSD[i,2]<-sd(matRMSE[,i])
}
  write.csv(matR_MeanSD, paste("GP/",model,"/Allpheno_GP_10-fold_gwSNP_MeanSD_RMSEvalue.csv", sep=""))

