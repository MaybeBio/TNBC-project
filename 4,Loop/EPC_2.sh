#!/bin/bash

#mamba activate bedpe
#首先进行简单的尝试，使用bed文件来注释bedPE文件
BedpeAnnotateFromBed -bed annotation.bed -in input.bedpe -out annotated_output.bedpe -col 5 -col_name GENE_NAME
BedpeAnnotateFromBed -bed 用于注释bedpe的bed文件like 启动子增强子文件等 \
                      -in 被注释的bedpe文件like loop文件 \
                      -out 输出的被注释文件 \
                      -col bed文件中哪一列用于注释，比如说是第4列，建议修改成SE1等等称呼  \
                      -col_name 注释列名字
                       
首先举个例子：
使用TNBC unique的loop数据，然后就是SE的数据
首先要对enhancer的文件进行一个标称,注意这里千万不能直接使用tail同名文件不然会变0
awk 'BEGIN {OFS="\t"} {print $0, "SE" NR-1}' TNBC_SE_sort.bed | tail -n+2 > TNBC_SE_sort_with_SE.bed

cut -f 1-6  loops_TNBC_unique_10000.bedpe  | tail -n+2 > loops_TNBC_unique_10000_16.bedpe

BedpeAnnotateFromBed -bed /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_SE_sort_with_SE.bed -in /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000.bedpe -out loops_TNBC_unique_10000_SEed.bedpe -col 6 -col_name SE_anno


如果还是使用bedtools intersect -a loop的bedPE文件（尝试使用anchor1） -b TE文件，SE文件，P文件（可以都输入） -wa -wb -filenames



