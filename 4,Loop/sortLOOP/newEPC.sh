#!/bin/bash

#�����ǻ�ȡunique��loop�����ݣ��������ݴ����Լ�id��ע
loopdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result
enhancerdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed
promoter=/mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19_promoter_sort.bed
myrep=(TNBC BT549 HCC70 MB231)


for rep in ${myrep[@]};do
  #�����ǻ�ȡunique��loop�����ݣ��������ݴ����Լ�id��ע
  awk 'BEGIN{OFS="\t"} {print $0, "loop"NR-1}' ${loopdir}/${rep}_preprocessing_any/loops_${rep}_unique_10000.bedpe > ${loopdir}/${rep}_preprocessing_any/loops_${rep}_unique_10000_withid.bedpe
  cut -f 1-3,12 ${loopdir}/${rep}_preprocessing_any/loops_${rep}_unique_10000_withid.bedpe | tail -n+2 > templeft.bed
  #����֮��ֻ����1-8�е������������Ԫ�����в��ÿ���
  bedtools intersect -a templeft.bed -b ${enhancerdir}/${rep}_SE_sort.bed ${enhancerdir}/${rep}_TE_sort.bed ${promoter} -wa -wb -loj -sortout -names SE TE P | cut -f 1-8 > ${rep}_${rep}_leftanchor.bed
rm templeft.bed
  #�����Ƕ�ͬһloop�е�HMEC����
  awk 'BEGIN{OFS="\t"} {print $0, "loop"NR-1}' ${loopdir}/${rep}_preprocessing_any/loops_HMEC_unique_10000.bedpe > ${loopdir}/${rep}_preprocessing_any/loops_HMEC_unique_10000_withid.bedpe
  cut -f 1-3,12 ${loopdir}/${rep}_preprocessing_any/loops_HMEC_unique_10000_withid.bedpe | tail -n+2 > templeft.bed
  bedtools intersect -a templeft.bed -b ${enhancerdir}/HMEC_SE_sort.bed ${enhancerdir}/HMEC_TE_sort.bed ${promoter} -wa -wb -loj -sortout -names SE TE P | cut -f 1-8 > ${rep}_HMEC_leftanchor.bed
rm templeft.bed

  #�����Ƕ���anchor����
  cut -f 4-6,12 ${loopdir}/${rep}_preprocessing_any/loops_${rep}_unique_10000_withid.bedpe | tail -n+2 > tempright.bed
  bedtools intersect -a tempright.bed -b ${enhancerdir}/${rep}_SE_sort.bed ${enhancerdir}/${rep}_TE_sort.bed ${promoter} -wa -wb -loj -sortout -names SE TE P | cut -f 1-8 > ${rep}_${rep}_rightanchor.bed
rm tempright.bed
  #�����Ƕ�ͬһloop�е�HMEC����
  cut -f 4-6,12 ${loopdir}/${rep}_preprocessing_any/loops_HMEC_unique_10000_withid.bedpe | tail -n+2 > tempright.bed
  bedtools intersect -a tempright.bed -b ${enhancerdir}/HMEC_SE_sort.bed ${enhancerdir}/HMEC_TE_sort.bed ${promoter} -wa -wb -loj -sortout -names SE TE P | cut -f 1-8 > ${rep}_HMEC_rightanchor.bed
rm tempright.bed
  
  python loop${rep}.py
done



#��һ�����ǿ�����κϲ�����anchor�е�ע�ͣ��Լ���κϲ�����anchor���˴�ѡ��anchor��ϵķ�ʽ���ϲ�anchor�����ݵ���anchor�ļ��е�loop id��
#��д��һ������