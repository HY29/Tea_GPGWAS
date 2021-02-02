#!/bin/bash
fastSTRUCTURE_DIR=${HOME}/miniconda3/pkgs/faststructure-1.0-py27h549429d_0/bin
PLINK_DIR=${HOME}/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/filtered_vcf
OUT_DIR=${HOME}/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/faststructure/list_gwas_maf_0.05_miss0.7

#Change Conda environment py27 (Python)
source activate py27
cd ${HOME}/miniconda3/pkgs/faststructure-1.0-py27h549429d_0/bin

python structure.py \
-K 3 \
--input=$PLINK_DIR/list_gwas_maf0.05_miss0.7 \
--output=$OUT_DIR/list_gwas_maf0.05_miss0.7 \
--format=bed 


