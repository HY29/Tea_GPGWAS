#!/bin/bash
PLINK_DIR=${HOME}/yamashita/tea_RAD-seq/tool
OUT_DIR=${HOME}/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/filtered_vcf

export PATH=$PLINK_DIR:$PATH

#plink command_make binary file (bed,bim,fam)
plink --file $OUT_DIR/list_gwas_maf0.05_miss0.7 \
--make-bed --out $OUT_DIR/list_gwas_maf0.05_miss0.7





