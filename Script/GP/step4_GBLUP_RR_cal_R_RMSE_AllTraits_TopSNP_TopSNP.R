#########################################################
### Calculation of R and RMSE value in GPwithGWAS
### Draft script 200811 HY
#########################################################

###Output of R & RMSE value in each round----
#set model
model = print(paste("GBLUP(RR)")) ##Caution!

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
  matR_Mean<-matrix(NA,nrow=length(nsnp))
  matR_SD<-matrix(NA,nrow=length(nsnp))
  for (n in 1:length(nsnp)){
    print(paste(z,n))
    one<-read.csv(paste("GPwithGWAS/",model,"/pheno_",z,"/Top_",nsnp[n],"/GP_10-fold_topSNP_Rvalue.csv", sep=""),row.names=1)
    matR[,n]<-one[,1]
    matR_Mean[n,]<-mean(matR[,n])
    matR_SD[n,] <-sd(matR[,n])
    write.csv(matR_Mean, paste("GPwithGWAS/",model,"/pheno_",z,"/AllTopSNP_GP_10-fold_topSNP_MeanRvalue.csv", sep=""))
    write.csv(matR_SD, paste("GPwithGWAS/",model,"/pheno_",z,"/AllTopSNP_GP_10-fold_topSNP_SDRvalue.csv", sep=""))
    }}

#Mean & SD of RMSE
for (z in 1:length(pheno)) {
  print(paste(n))
  matRMSE<-matrix(NA,nrow=length(R.value),ncol=length(nsnp))
  matRMSE_Mean<-matrix(NA,nrow=length(nsnp))
  matRMSE_SD<-matrix(NA,nrow=length(nsnp))
  for (n in 1:length(nsnp)){
    print(paste(z,n))
    one<-read.csv(paste("GPwithGWAS/",model,"/pheno_",z,"/Top_",nsnp[n],"/GP_10-fold_topSNP_RMSEvalue.csv", sep=""),row.names=1)
    matRMSE[,n]<-one[,1]
    matRMSE_Mean[n,]<-mean(matRMSE[,n])
    matRMSE_SD[n,] <-sd(matRMSE[,n])
    write.csv(matR_Mean, paste("GPwithGWAS/",model,"/pheno_",z,"/AllTopSNP_GP_10-fold_topSNP_MeanRMSEvalue.csv", sep=""))
    write.csv(matR_Mean, paste("GPwithGWAS/",model,"/pheno_",z,"/AllTopSNP_GP_10-fold_topSNP_SDRMSEvalue.csv", sep=""))
    }}
