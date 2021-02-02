#!/bin/bash

cd ${HOME}/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/filtered_vcf

ScriptDIR=${HOME}/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/Script

#Rscript
Rscript --vanilla $ScriptDIR/SNPImp_missForest.R

