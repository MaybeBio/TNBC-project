#!/bin/bash
# Script to compute CTCF orientation using Maggie and Pratiik's CTCF orientation tool at https://github.com/magmarshh/CTCF_orientation 

# Files, same for each data type 
#FILE=( HMEC_unique BT549_unique HMEC_common BT549_common )
# hg19 fimo file,见上使用的是peak calling结果文件
#FIMODIR=/mnt/disk4/haitao/bysj_seu/geo_data/other_omics/CTCF 

# Mustache Loops Data Directory where the loop files are stored 
#LOOPDIR=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/preprocessing_any
# Mustache results directory 
#DIROUT=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/preprocessing_any/Loops
# Loop resolution
#LOOPRES=( 10kb )

# Loop through every condition and resolution and produce the CTCF piechart
#for f in ${FILE[@]}; do
#		python /mnt/disk4/haitao/software/CTCF_orientation/ctcf_orientation.py -l ${LOOPDIR}/loops_${f}_10000.bedpe \
#		-m ${FIMODIR}/${f%%_*}*.bed -o ${DIROUT}/Mustache_${f}_10000_CTCF_orientation.pdf
#done

#问题就是这里如何使用对应的bed文件，因为对应的loop文件是不一致的，但是所提供的motif bed文件是使用peak calling出来的结果，不是统一公用的hg19的motif文件，所以如何对应上
#要么就是分开来分析
#因为每次使用的file是HMEC_unique的，所以可以直接去掉后缀取HMEC，然后注意到CTCF文件夹中的motif文件名是HMEC_rep2_summits.bed，如果可以模糊匹配的话，就使用file-去后缀-*-.bed
#${file%%_*}就是

#上面都是旧代码
# hg19 CTCF fimo.bed file will be used 注意这里使用的也是对应rep的chip-seq peak calling的结果peak文件，没有使用公共motif分析预测的bed文件——
#1，首先使用的肯定是合并之后的每个sample的motif的bed文件，
#2，然后macs2其实是能够使用多rep作为input的，即可以做到TNBC 3合1的
#暂时使用这个"/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/hg19_CTCF.bed"

myrep=(TNBC) #HMEC BT549 HCC70 MB231
CTCFdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/CTCFbed
CTCFbed=(CTCF CTCF_hg19_all_filter_reduce CTCF_hg19_all)      #bed文件也可以列一个数组然后执行  
loopdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result  #loop
# File suffix, corresponds to the loop conditions 
  declare -A suffix=(
    ["TNBC"]="TNBC_unique TNBC_common HMEC_unique HMEC_common TNBC_all HMEC_all"
)
#下面是原来的common+unique,现在补充上面的all
#["BT549"]="BT549_unique BT549_common HMEC_unique HMEC_common BT549_all HMEC_all"
#["HCC70"]="HCC70_unique HCC70_common HMEC_unique HMEC_common HCC70_all HMEC_all"
#["MB231"]="MB231_unique MB231_common HMEC_unique HMEC_common MB231_all HMEC_all"

for bed in ${CTCFbed[@]};do
  for rep in ${myrep[@]};do
  fs_array=(${suffix[$rep]})
   for fs in  ${fs_array[@]};do
    python /mnt/disk4/haitao/software/CTCF_orientation/ctcf_orientation.py -l ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe -m ${CTCFdir}/${bed}.bed -o ${loopdir}/${rep}_preprocessing_any/${bed}_${fs}_10kb_CTCF_orientation.pdf
   done;
 done;
done


