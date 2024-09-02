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
  #先处理HMEC的循环，将识别出来的文件
#  java -Xms10G -Xmx20G -jar /mnt/disk4/haitao/software/juicer_tools.2.20.00.jar motifs hg19 ${bed_file_dir}/HMEC ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe
#  done;
#  mkdir -p /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/${rep}
#  mv *with_motifs* /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/${rep}
  
  #其次是使用自身的peak文件
  for fs in  ${fs_array[@]};do
  java -Xms5G -Xmx10G -jar /mnt/disk4/haitao/software/juicer_tools.2.18.00.jar motifs hg19 ${bed_file_dir}/${rep} ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe
  done;
done


#主要是命令还有一个选项，就是使用motif文件可以指定获取motif文件，这里可以和前面CTCF方向计算中的bed文件相统一
#这里首先不清楚各自TNBC中是否有必要再使用HMEC的bed文件；其次，是如何使用HMEC的peak文件去注释，指定输出不重复的话;
#注意各个文件的甄别   
#HMEC先不用