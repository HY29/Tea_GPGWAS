#!/bin/bash

ToolDir=$HOME/yamashita/tea_RAD-seq/tool
POPDir=$HOME/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/popmap
STACKS_Dir=$HOME/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/stacks
OUT_Dir=$HOME/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/vcf


export PATH=$ToolDir/stacks-2.5/bin:$PATH


#stacks_populations command_pipeline
populations -P $STACKS_Dir/ -M $POPDir/rad_all_popmap.txt -O $OUT_Dir/ -t 8 --vcf 


