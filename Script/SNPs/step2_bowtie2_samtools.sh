#!/bin/bash
RefDir=/home/teaplant/TeaGenome_ChrLev
DataDir=/home/teaplant/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/samples
Bowtie2_bin=$HOME/yamashita/tea_RAD-seq/tool/bowtie2-2.3.5.1-linux-x86_64
Samtools_bin=/usr/local/samtools-1.9

export PATH=$Bowtie2_bin:$Samtools_bin:$PATH

mkdir ../sortBAM
cd ../sortBAM

#CsiRAD3
#Mapping SE
for i in $(seq -w 1 184)
do

#bowtie2_samtools
bowtie2 -x $RefDir/CSS_ChrLev_20200506_Genome \
-U $DataDir/sample${i}.fq \
-S sample${i}.sam

samtools sort \
  -@ 8 \
  -o sample${i}.bam \
  sample${i}.sam
  
  rm sample${i}.sam
  
done
