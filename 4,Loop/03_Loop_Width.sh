#!/bin/bash
# Script to use LoopWidth https://github.com/magmarshh/LoopWidth to compute loop widths for 

# Conditions, correspond to file names 
#FILE1=( HMEC_unique BT549_unique HMEC_common BT549_common HMEC_all BT549_all) #BT549
#FILE2=( HMEC_unique HCC70_unique HMEC_common HCC70_common HMEC_all HCC70_all) #HCC70
#FILE3=( HMEC_unique MB231_unique HMEC_common MB231_common HMEC_all MB231_all) #MB231
FILE=(HMEC_unique TNBC_unique HMEC_common TNBC_common HMEC_all TNBC_all) #TNBC

# Mustache Loops Data Directory，用于存放00脚本中执行结果的bed，bedpe文件的地方
#LOOPDIR1=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/BT549_preprocessing_any #BT549
#LOOPDIR2=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/HCC70_preprocessing_any #HCC70
#LOOPDIR3=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/MB231_preprocessing_any #MB231
LOOPDIR=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any #TNBC

# Results directory 结果文件
#DIROUT=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/preprocessing_any/Loops
# Loop resolution，可以使用20kb以及40kb的数据都计算一次，下面主要是看loop实际的分辨率有多少个，才会使用第一层循环，即res in 1-4，但是现在只有1个10kb的分辨率，其实可以不设置这一层循环
#LOOPRES=(10kb)
#LOOPNUM=(10000)

# Loop through every condition and resolution and output loop piecharts，res之后再使用{1..4}，下面的代码都进行了修改，只有一个文件的index都尽量不使用了
#因为只有一个res参数，所以直接使用[1]，后续使用多分辨的时候，就在最外面再套一层loop：
#for res in {1..2};do
#对应 ${LOOPNUM[$res]}
#done;

	for f in ${FILE[@]}; do
		python /mnt/disk4/haitao/software/LoopWidth/LoopWidth_piechart.py --loop ${LOOPDIR}/loops_${f}_10000.bedpe --res 10000 \
		--output ${LOOPDIR}/Mustache_${f}_10kb_loop_sizes_piechart.pdf;
   python /mnt/disk4/haitao/software/LoopWidth/LoopWidth_violinplot.py --loop ${LOOPDIR}/loops_${f}_10000.bedpe --labels ${f}_10000  \
		--output ${LOOPDIR}/Mustache_${f}_10kb_loop_sizes_violinplot.pdf;
	done
  
#	for f in ${FILE1[@]}; do
#		python /mnt/disk4/haitao/software/LoopWidth/LoopWidth_piechart.py --loop ${LOOPDIR1}/loops_${f}_10000.bedpe --res 10000 \
#		--output ${LOOPDIR1}/Mustache_${f}_10kb_loop_sizes_piechart.pdf;
#   python /mnt/disk4/haitao/software/LoopWidth/LoopWidth_violinplot.py --loop ${LOOPDIR1}/loops_${f}_10000.bedpe --labels ${f}_10000  \
#		--output ${LOOPDIR1}/Mustache_${f}_10kb_loop_sizes_violinplot.pdf;
#	done
	
#  for f in ${FILE2[@]}; do
#		python /mnt/disk4/haitao/software/LoopWidth/LoopWidth_piechart.py --loop ${LOOPDIR2}/loops_${f}_10000.bedpe --res 10000 \
#		--output ${LOOPDIR2}/Mustache_${f}_10kb_loop_sizes_piechart.pdf;
#   python /mnt/disk4/haitao/software/LoopWidth/LoopWidth_violinplot.py --loop ${LOOPDIR2}/loops_${f}_10000.bedpe --labels ${f}_10000  \
#		--output ${LOOPDIR2}/Mustache_${f}_10kb_loop_sizes_violinplot.pdf;
#	done
 
# 	for f in ${FILE3[@]}; do
#		python /mnt/disk4/haitao/software/LoopWidth/LoopWidth_piechart.py --loop ${LOOPDIR3}/loops_${f}_10000.bedpe --res 10000 \
#		--output ${LOOPDIR3}/Mustache_${f}_10kb_loop_sizes_piechart.pdf;
#   python /mnt/disk4/haitao/software/LoopWidth/LoopWidth_violinplot.py --loop ${LOOPDIR3}/loops_${f}_10000.bedpe --labels ${f}_10000  \
#		--output ${LOOPDIR3}/Mustache_${f}_10kb_loop_sizes_violinplot.pdf;
#	done


