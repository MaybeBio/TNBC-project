#!/bin/bash

#首先是获取unique的loop的数据，做好数据处理以及id标注
loopdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result
enhancerdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed
promoter=/mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19_promoter_sort.bed
myrep=(TNBC BT549 HCC70 MB231)


for rep in ${myrep[@]};do
  #首先是获取unique的loop的数据，做好数据处理以及id标注
  awk 'BEGIN{OFS="\t"} {print $0, "loop"NR-1}' ${loopdir}/${rep}_preprocessing_any/loops_${rep}_unique_10000.bedpe > ${loopdir}/${rep}_preprocessing_any/loops_${rep}_unique_10000_withid.bedpe
  cut -f 1-3,12 ${loopdir}/${rep}_preprocessing_any/loops_${rep}_unique_10000_withid.bedpe | tail -n+2 > templeft.bed
  #交集之后只保留1-8列的数据其余调控元件的列不用看了
  bedtools intersect -a templeft.bed -b ${enhancerdir}/${rep}_SE_sort.bed ${enhancerdir}/${rep}_TE_sort.bed ${promoter} -wa -wb -loj -sortout -names SE TE P | cut -f 1-8 > ${rep}_${rep}_leftanchor.bed
rm templeft.bed
  #下面是对同一loop中的HMEC而言
  awk 'BEGIN{OFS="\t"} {print $0, "loop"NR-1}' ${loopdir}/${rep}_preprocessing_any/loops_HMEC_unique_10000.bedpe > ${loopdir}/${rep}_preprocessing_any/loops_HMEC_unique_10000_withid.bedpe
  cut -f 1-3,12 ${loopdir}/${rep}_preprocessing_any/loops_HMEC_unique_10000_withid.bedpe | tail -n+2 > templeft.bed
  bedtools intersect -a templeft.bed -b ${enhancerdir}/HMEC_SE_sort.bed ${enhancerdir}/HMEC_TE_sort.bed ${promoter} -wa -wb -loj -sortout -names SE TE P | cut -f 1-8 > ${rep}_HMEC_leftanchor.bed
rm templeft.bed

  #下面是对右anchor而言
  cut -f 4-6,12 ${loopdir}/${rep}_preprocessing_any/loops_${rep}_unique_10000_withid.bedpe | tail -n+2 > tempright.bed
  bedtools intersect -a tempright.bed -b ${enhancerdir}/${rep}_SE_sort.bed ${enhancerdir}/${rep}_TE_sort.bed ${promoter} -wa -wb -loj -sortout -names SE TE P | cut -f 1-8 > ${rep}_${rep}_rightanchor.bed
rm tempright.bed
  #下面是对同一loop中的HMEC而言
  cut -f 4-6,12 ${loopdir}/${rep}_preprocessing_any/loops_HMEC_unique_10000_withid.bedpe | tail -n+2 > tempright.bed
  bedtools intersect -a tempright.bed -b ${enhancerdir}/HMEC_SE_sort.bed ${enhancerdir}/HMEC_TE_sort.bed ${promoter} -wa -wb -loj -sortout -names SE TE P | cut -f 1-8 > ${rep}_HMEC_rightanchor.bed
rm tempright.bed
  
  python loop${rep}.py
done



#下一步就是考虑如何合并单个anchor中的注释，以及如何合并两个anchor，此处选用anchor组合的方式来合并anchor，依据的是anchor文件中的loop id列
#编写了一个程序