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



#对于TNBC，可以合并来处理，只处理unique的loop文件
cut -f 1-6  loops_TNBC_unique_10000.bedpe  | tail -n+2 > loops_TNBC_unique_10000_16.bedpe
bedtools intersect -a /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000_16.bedpe  -b /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_SE_sort.bed  /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_TE_sort.bed  /mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19_promoter_sort.bed -wa -wb -loj -sortout -names SE TE P
                       
                       
                       
                       

首先是在loop文件中增加一个标识符
awk 'BEGIN{OFS="\t"} {print $0, "loop"NR-1}' /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000.bedpe > loops_TNBC_unique_10000_withid.bedpe
或者使用临时文件中转，不用次次转换:至少对于每一个bed文件，都只取对应的bed3格式列，然后一般第一行都去掉
awk 'BEGIN{OFS="\t"} {print $0, "loop"NR-1}' /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000.bedpe >  /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000_withid.bedpe
cut -f 1-3,12 /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000_withid.bedpe | tail -n+2 > temp.bed
bedtools intersect -a temp.bed -b /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_SE_sort.bed /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_TE_sort.bed /mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19_promoter_sort.bed -wa -wb -loj -sortout -names SE TE P
rm temp.bed

cut -f 4-6,12 /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000_withid.bedpe | tail -n+2 > temp.bed
bedtools intersect -a temp.bed -b /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_SE_sort.bed /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_TE_sort.bed /mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19_promoter_sort.bed -wa -wb -loj -sortout -names SE TE P
rm temp.bed

现在有一个问题，就是如何将两者的数据和在一起，必须按照loop配对的顺序

#或者是3者分开来处理，获取都被SE调控的gene
也可以对HMEC做一个分析，最后将TNBC中SE的EPC-HMEC的unique中的SE-PC的gene作为最后模型输入



