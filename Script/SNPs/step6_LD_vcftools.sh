#!/bin/bash
VCFTOOLS_DIR=${HOME}/yamashita/tea_RAD-seq/tool/vcftools-v0.1.16/bin
Filltered_VCF_DIR=/home/teaplant/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/filtered_vcf
LD_DIR=/home/teaplant/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/LD

export PATH=$VCFTOOLS_DIR:$PATH

#vcftools command
vcftools \
--vcf $Filltered_VCF_DIR/list_gwas_maf0.05_miss0.7.recode.vcf \
--geno-r2 \
--ld-window-bp 50000 \
--out $LD_DIR/list_gwas_maf0.05_miss0.7_ld_window_50kb




