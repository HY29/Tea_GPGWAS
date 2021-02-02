#!/bin/bash
STACKS_DIR=${HOME}/yamashita/tea_RAD-seq/tool/stacks-2.5/bin
ToolDir=$HOME/yamashita/tea_RAD-seq/tool
POPDir=$HOME/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/popmap
BAMDir=$HOME/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/sortBAM
OUTDir=$HOME/yamashita/tea_RAD-seq/2018_tea_accessions_shizuoka/CsiRAD3/ChrLev/stacks

export PATH=$ToolDir/stacks-2.5/bin:$PATH

cd ../stacks

#stacks_pipeline
ref_map.pl -T 4 \
--samples $BAMDir/ \
--popmap $POPDir/rad_all_popmap.txt \
-o $OUTDir/ \
