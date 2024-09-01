#!/bin/bash

#mamba activate HiCExplorer
myrep=(HMEC_rep1 HMEC_rep2 BT549_rep1 BT549_rep2 HCC70_rep1 HCC70_rep2 MB231)
folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script4

for rep in ${myrep[@]};do

  hicMergeDomains --domainFiles ${folder}/10kb/10000_${rep}_domains.bed ${folder}/20kb/20000_${rep}_domains.bed ${folder}/40kb/40000_${rep}_domains.bed ${folder}/100kb/100000_${rep}_domains.bed ${folder}/500kb/500000_${rep}_domains.bed  --outputMergedList  ${rep}_merged_noCTCF --outputRelationList  ${rep}_relation_noCTCF  --outputTreePlotPrefix ${rep}_plot_noCTCF  --outputTreePlotFormat pdf

done

#�CproteinFile ctcf_sorted.bed Ϊ���ܹ����õ�����TADs֮��Ĺ�ϵ���ɰ�����صĵ������ļ�(�������ڲ��鶯���CTCF),��Ҫʹ��broadpeak��ʽ��peak�ļ�������ÿһ��rep����peak�����ʹ��idr֮��һ���Ե�peak���������������TAD֮��Ĺ�ϵ�����ʹ�õ���rep��1peak�͸������ˣ�����ʹ��һ���Եĵ�peak֮�⣬������ʹ�ù������ݿ��е�bed�ļ�������˵�ǹ�����CTCF��hg19��bed�ļ�




 
 
 
 