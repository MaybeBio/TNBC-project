#!/bin/bash

#mamba activate HiCExplorer
myrep=(HMEC_rep1 HMEC_rep2 BT549_rep1 BT549_rep2 HCC70_rep1 HCC70_rep2 MB231)
folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script4

for rep in ${myrep[@]};do

  hicMergeDomains --domainFiles ${folder}/10kb/10000_${rep}_domains.bed ${folder}/20kb/20000_${rep}_domains.bed ${folder}/40kb/40000_${rep}_domains.bed ${folder}/100kb/100000_${rep}_domains.bed ${folder}/500kb/500000_${rep}_domains.bed  --outputMergedList  ${rep}_merged_noCTCF --outputRelationList  ${rep}_relation_noCTCF  --outputTreePlotPrefix ${rep}_plot_noCTCF  --outputTreePlotFormat pdf

done

#–proteinFile ctcf_sorted.bed 为了能够更好地评估TADs之间的关系，可包括相关的蛋白质文件(例如用于哺乳动物的CTCF),需要使用broadpeak格式的peak文件，但是每一个rep都有peak，如果使用idr之后一致性的peak，那如何评估各个TAD之间的关系，如果使用单个rep的1peak就更不行了；除了使用一致性的的peak之外，还可以使用公共数据库中的bed文件，比如说是公共的CTCF的hg19的bed文件




 
 
 
 