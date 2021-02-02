#!/bin/bash
ToolDir=$HOME/yamashita/tea_RAD-seq/tool

export PATH=$ToolDir/bowtie2-2.3.5.1-linux-x86_64:$PATH

cd /home/teaplant/TeaGenome_ChrLev

bowtie2-build -f CSS_ChrLev_20200506_Genome.fas CSS_ChrLev_20200506_Genome

