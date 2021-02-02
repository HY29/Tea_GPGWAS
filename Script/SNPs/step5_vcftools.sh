#!/bin/bash
VCFTOOLS_DIR=${HOME}/yamashita/tea_RAD-seq/tool/vcftools-v0.1.16/bin
VCF_DIR=/home/teaplant/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/vcf
LIST_DIR=/home/teaplant/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/samplelist_vcffilter
OUT_DIR=/home/teaplant/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/filtered_vcf

export PATH=$VCFTOOLS_DIR:$PATH

#vcftools command
vcftools \
--vcf $VCF_DIR/populations.snps.vcf \
--keep $LIST_DIR/list_gwas.txt \
--recode --out $OUT_DIR/list_gwas_maf0.05_miss0.7 \
--maf 0.05 \
--max-missing 0.7 \




