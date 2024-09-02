#!/bin/bash

myrep=(BT549 HCC70 MB231) #HMEC BT549 HCC70 MB231
hicdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/hic2juice/merge_fq #hic
loopdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result  #loop
# File suffix, corresponds to the loop conditions 
HMECsuffix=( HMEC_unique HMEC_common )
BT549suffix=( BT549_unique BT549_common )
HCC70suffix=( HCC70_unique HCC70_common )
MB231suffix=( MB231_unique MB231_common )

#���淽���в�ͨ�����������Զ������飺
declare -A suffix=(
    ["BT549"]="BT549_unique BT549_common BT549_all HMEC_unique HMEC_common HMEC_all"
    ["HCC70"]="HCC70_unique HCC70_common HCC70_all HMEC_unique HMEC_common HMEC_all" 
    ["MB231"]="MB231_unique MB231_common MB231_all HMEC_unique HMEC_common HMEC_all"
)
  
for rep in ${myrep[@]};do
  fs_array=(${suffix[$rep]})
  for fs in  ${fs_array[@]};do
    java -Xms10G -Xmx20G -jar /mnt/disk4/haitao/software/juicer_tools.2.20.00.jar apa -k NONE -u -r 10000 \
  ${hicdir}/${rep}/${rep}.allValidPairs.hic  ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe \
  ${loopdir}/${rep}_preprocessing_any/perchrAPA/${rep}_${fs} 

    java -Xms10G -Xmx20G -jar /mnt/disk4/haitao/software/juicer_tools.2.20.00.jar apa -k NONE  -u -r 10000 \
  ${hicdir}/HMEC/HMEC.allValidPairs.hic  ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe \
  ${loopdir}/${rep}_preprocessing_any/perchrAPA/HMEC_${fs}
  done;
done
  
  
#loop��peak�ļ��ͷֿ���bedpe�ļ������ֿ���unique��common�ļ���ʲô�����𣬿��Ժϲ���һ�������
#ǰ���hic��loop�ļ�Ӧ���Ƕ�Ӧ��
#����ļ���${rep}_${fs}��ǰ����hic��������loop�ļ�
#���������ԭʼ����-u�ģ�������-u֮��call�����ľ���perchr��������ȫ��

#������δ����ű�
#for rep in ${myrep[@]};do
#  fs_array=(${suffix[$rep]})
#  for fs in  ${fs_array[@]};do
#    java -Xms5G -Xmx10G -jar /mnt/disk4/haitao/software/juicer_tools.2.20.00.jar apa -k NONE -r 10000 \
#  ${hicdir}/${rep}/${rep}.allValidPairs.hic  ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe \
#  ${loopdir}/${rep}_preprocessing_any/${fs}_new 
#  done;
#  for fs in ${HMECsuffix[@]};do
#    java -Xms5G -Xmx10G -jar /mnt/disk4/haitao/software/juicer_tools.2.20.00.jar apa -k NONE  -r 10000 \
#  ${hicdir}/HMEC/HMEC.allValidPairs.hic  ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe \
#  ${loopdir}/${rep}_preprocessing_any/${fs}_new
#  done;
#done