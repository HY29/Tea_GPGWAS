###vcf data construction for GWAS

library(vcfR)

#read_compressed_vcffile_*.vcf.gz
vcf <- read.vcfR("../filtered_vcf/list_gwas_maf0.05_miss0.7.recode.vcf")
#extract genotype data from vcf
gt <- extract.gt(x=vcf, element = "GT", IDtoRowNames = FALSE)
#get marker information
chr <- getCHROM(vcf)
pos <- getPOS(vcf)
id <- getID(vcf)

#create a matrix of gt score
gt.score <- matrix(NA, nrow(gt), ncol(gt))
gt.score[gt == "0/0"] <- -1
gt.score[gt == "0/1"] <- 0
gt.score[gt == "1/1"] <- 1

#name the rows and columns of matrix
rownames(gt.score) <- rownames(gt)
colnames(gt.score) <- colnames(gt)

#transpose the matrix
gt.score <- t(gt.score)


#format to dataframe
chr <- as.data.frame(chr)
pos <- as.data.frame(pos)
ID <- as.data.frame(id)
gt.score <- as.data.frame(gt.score)
gt.score <- t(gt.score)
SNPtable <- cbind.data.frame(chr, pos, ID, gt.score)
snp_only <- SNPtable[4:153]

###snp imputation by missForest
library(doParallel)
library(iterators)
library(missForest)
registerDoParallel(cores=8)

snp.imp <- missForest(snp_only, verbose=TRUE, parallelize="variables")
snp.imp.ximp <- snp.imp$ximp
snp.table <- cbind.data.frame(chr, pos, ID, snp.imp.ximp)

write.csv(snp.table, "../filtered_vcf/list_gwas_maf0.05_miss0.7.imp.csv")
