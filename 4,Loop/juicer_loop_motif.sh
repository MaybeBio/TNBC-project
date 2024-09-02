#!/bin/bash

myrep=(BT549 HCC70 MB231) #HMEC BT549 HCC70 MB231
loopdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result  #loop
bed_file_dir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/CTCFbed/peak_bed_track
declare -A suffix=(
    ["BT549"]="BT549_unique BT549_common BT549_all HMEC_unique HMEC_common HMEC_all"
    ["HCC70"]="HCC70_unique HCC70_common HCC70_all HMEC_unique HMEC_common HMEC_all" 
    ["MB231"]="MB231_unique MB231_common MB231_all HMEC_unique HMEC_common HMEC_all" 
)

for rep in ${myrep[@]};do
  fs_array=(${suffix[$rep]})
  
#  for fs in  ${fs_array[@]};do
  #�ȴ���HMEC��ѭ������ʶ��������ļ�
#  java -Xms10G -Xmx20G -jar /mnt/disk4/haitao/software/juicer_tools.2.20.00.jar motifs hg19 ${bed_file_dir}/HMEC ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe
#  done;
#  mkdir -p /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/${rep}
#  mv *with_motifs* /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/${rep}
  
  #�����ʹ�������peak�ļ�
  for fs in  ${fs_array[@]};do
  java -Xms5G -Xmx10G -jar /mnt/disk4/haitao/software/juicer_tools.2.18.00.jar motifs hg19 ${bed_file_dir}/${rep} ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe
  done;
done


#��Ҫ�������һ��ѡ�����ʹ��motif�ļ�����ָ����ȡmotif�ļ���������Ժ�ǰ��CTCF��������е�bed�ļ���ͳһ
#�������Ȳ��������TNBC���Ƿ��б�Ҫ��ʹ��HMEC��bed�ļ�����Σ������ʹ��HMEC��peak�ļ�ȥע�ͣ�ָ��������ظ��Ļ�;
#ע������ļ������   
#HMEC�Ȳ���