#!/bin/bash

loopdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result
myrep=(BT549 HCC70 MB231 TNBC) #HMEC BT549 HCC70 MB231
declare -A suffix=(
    ["BT549"]="BT549_unique BT549_common BT549_all HMEC_unique HMEC_common HMEC_all"
    ["HCC70"]="HCC70_unique HCC70_common HCC70_all HMEC_unique HMEC_common HMEC_all" 
    ["MB231"]="MB231_unique MB231_common MB231_all HMEC_unique HMEC_common HMEC_all" 
    ["TNBC"]="TNBC_unique TNBC_common TNBC_all HMEC_unique HMEC_common HMEC_all"
)


for rep in ${myrep[@]};do
  fs_array=(${suffix[$rep]})
  for fs in  ${fs_array[@]};do
    
# 提取第一列、第二列和第五列，并将第二列和第五列加上5000，并修改列名
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print "chr", "anchor1", "anchor2"} NR>1 {print $1, $2+5000, $5+5000}' ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe > ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000_summit.bed

  done;
done