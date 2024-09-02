#!/bin/bash
# Script to compute CTCF orientation using Maggie and Pratiik's CTCF orientation tool at https://github.com/magmarshh/CTCF_orientation 

# Files, same for each data type 
#FILE=( HMEC_unique BT549_unique HMEC_common BT549_common )
# hg19 fimo file,����ʹ�õ���peak calling����ļ�
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

#��������������ʹ�ö�Ӧ��bed�ļ�����Ϊ��Ӧ��loop�ļ��ǲ�һ�µģ��������ṩ��motif bed�ļ���ʹ��peak calling�����Ľ��������ͳһ���õ�hg19��motif�ļ���������ζ�Ӧ��
#Ҫô���Ƿֿ�������
#��Ϊÿ��ʹ�õ�file��HMEC_unique�ģ����Կ���ֱ��ȥ����׺ȡHMEC��Ȼ��ע�⵽CTCF�ļ����е�motif�ļ�����HMEC_rep2_summits.bed���������ģ��ƥ��Ļ�����ʹ��file-ȥ��׺-*-.bed
#${file%%_*}����

#���涼�Ǿɴ���
# hg19 CTCF fimo.bed file will be used ע������ʹ�õ�Ҳ�Ƕ�Ӧrep��chip-seq peak calling�Ľ��peak�ļ���û��ʹ�ù���motif����Ԥ���bed�ļ�����
#1������ʹ�õĿ϶��Ǻϲ�֮���ÿ��sample��motif��bed�ļ���
#2��Ȼ��macs2��ʵ���ܹ�ʹ�ö�rep��Ϊinput�ģ�����������TNBC 3��1��
#��ʱʹ�����"/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/hg19_CTCF.bed"

myrep=(TNBC) #HMEC BT549 HCC70 MB231
CTCFdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/CTCFbed
CTCFbed=(CTCF CTCF_hg19_all_filter_reduce CTCF_hg19_all)      #bed�ļ�Ҳ������һ������Ȼ��ִ��  
loopdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result  #loop
# File suffix, corresponds to the loop conditions 
  declare -A suffix=(
    ["TNBC"]="TNBC_unique TNBC_common HMEC_unique HMEC_common TNBC_all HMEC_all"
)
#������ԭ����common+unique,���ڲ��������all
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


