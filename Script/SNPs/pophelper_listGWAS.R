#########################################################
### Barplot of fastSTRUCTURE results
### Draft script 200313 HY
#########################################################

#load library
library(pophelper)
library(ggplot2)
library(ggpubr)
library(gridExtra)

#setwd
setwd("/Users/HirotoYAMASHITA/BioInfo/Ranalysis/Genetics/TeaGWAS/ChrLev")

#check version
packageDescription("pophelper", fields="Version")
#input more run files
sfiles <- list.files("faststructure/")
#input individuals information
inds <- read.csv("150sample_group.csv", header = TRUE,stringsAsFactors = F)
#input onelabel
indsinf <- inds[,2,drop=FALSE]
#convert q-matrix run files
setwd("faststructure")
slist <- readQ(files=sfiles, 
               filetype = "auto"
               )
#if all runs are equal length, add indlab to all runs
if(length(unique(sapply(slist,nrow)))==1) slist <- lapply(slist,"rownames<-",inds$sample)

#setting of font
library(showtext)
font_add_google("Patrick Hand", "patrick")
showtext_auto()

#list of color variation
clist <- list(
  "shiny"=c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E"),
  "strong"=c("#11A4C8","#63C2C5","#1D4F9F","#0C516D","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#ED3995","#7E277C","#F7EC16","#F8941E","#8C2A1C","#808080"),
  "oceanfive"=c("#00A0B0", "#6A4A3C", "#CC333F", "#EB6841", "#EDC951"),
  "keeled"=c("#48B098", "#91CB62", "#FFEE3B", "#FB9013", "#FF3C28"),
  "vintage"=c("#400F13", "#027368", "#A3BF3F", "#F2B950", "#D93A2B"),
  "muted"=c("#46BDDD","#82DDCE","#F5F06A","#F5CC6A","#F57E6A"),
  "teal"=c("#CFF09E","#A8DBA8","#79BD9A","#3B8686","#0B486B"),
  "merry"=c("#5BC0EB","#FDE74C","#9BC53D","#E55934","#FA7921"),
  "funky"=c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762","#ED8F47","#9471B4"),
  "retro"=c("#01948E","#A9C4E2","#E23560","#01A7B3","#FDA963","#323665","#EC687D"),
  "cb_paired"=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
  "cb_set3"=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
  "morris"=c("#4D94CC","#34648A","#8B658A","#9ACD32","#CC95CC","#9ACD32","#8B3A39","#CD6601","#CC5C5B","#8A4500"),
  "wong"=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#006699","#D55E00","#CC79A7"),
  "krzywinski"=c("#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE")
  )

#slist[c(12,14:20,2:6)] #K=2-14

#ploting
p <- plotQ(slist[c(1:2)],
           returnplot = T,
           exportplot = T,
           quiet = T,
           showyaxis = T,
           showindlab = F,
           indlabspacer = 0.1,
           clustercol = clist$shiny, # Setting of color
           sortind = "label",
           grplab = indsinf,
           ordergrp = T,
           showlegend = F, #Setting of legend TRUE or FALSE
           legendkeysize = 3,
           legendtextsize = 8,
           legendpos = "right",
           #selgrp = "site", # setting on the label based on plotting
           #grplabpos = 0, #default = 0
           #subset=c("Kyoto"),
           #font="patrick",
           splab = paste0("K=", sapply(slist[c(1:2)],ncol)),
           splabsize = 10,
           splabcol = "black",
           linesize = 0.3,
           linecol = "black",
           pointsize = 3,
           pointcol = "black",
           grplabsize = 3,
           grplabangle = 0,
           grplabcol = "black",
           grplabheight = 1,
           grplabspacer = 0.01,
           panelspacer = 0.3,
           basesize=12,
           imgoutput = "join",
           height = 1,
           width = 8,
           dpi = 350,
           imgtype = "tiff",
           #outputfilename = "Demo",
           )





