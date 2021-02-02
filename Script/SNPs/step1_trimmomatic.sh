#!/bin/bash

DataDir=$HOME/yamashita/tea_RAD-seq/2019_tea_accessions_Japan/raw
ToolDir=$HOME/yamashita/tea_RAD-seq/tool
Trimmomatic_jar=$ToolDir/Trimmomatic-0.39/trimmomatic-0.39.jar
Trimmomatic_adapters=$ToolDir/Trimmomatic-0.39/adapters

mkdir ../after_trimmomatic
cd ../after_trimmomatic

#CsiRAD4
for i in $(seq -w 1 184)
do
# Trimmomatic
  java -jar $Trimmomatic_jar SE \
    -threads 4 -phred33 -trimlog Trimmomatic_CsiRAD3_d1${i}.log \
    $DataDir/CsiRAD3_d1${i}.fastq.gz \
    sample0${i}.fastq.gz \
    ILLUMINACLIP:$Trimmomatic_adapters/TruSeq3-PE-2.fa:2:30:10 \
    LEADING:19 TRAILING:19 SLIDINGWINDOW:30:20 AVGQUAL:20 MINLEN:51
done



