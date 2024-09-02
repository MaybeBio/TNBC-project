#!/bin/bash

#mamba activate bedpe
#���Ƚ��м򵥵ĳ��ԣ�ʹ��bed�ļ���ע��bedPE�ļ�
BedpeAnnotateFromBed -bed annotation.bed -in input.bedpe -out annotated_output.bedpe -col 5 -col_name GENE_NAME
BedpeAnnotateFromBed -bed ����ע��bedpe��bed�ļ�like ��������ǿ���ļ��� \
                      -in ��ע�͵�bedpe�ļ�like loop�ļ� \
                      -out ����ı�ע���ļ� \
                      -col bed�ļ�����һ������ע�ͣ�����˵�ǵ�4�У������޸ĳ�SE1�ȵȳƺ�  \
                      -col_name ע��������
                       
���Ⱦٸ����ӣ�
ʹ��TNBC unique��loop���ݣ�Ȼ�����SE������
����Ҫ��enhancer���ļ�����һ�����,ע������ǧ����ֱ��ʹ��tailͬ���ļ���Ȼ���0
awk 'BEGIN {OFS="\t"} {print $0, "SE" NR-1}' TNBC_SE_sort.bed | tail -n+2 > TNBC_SE_sort_with_SE.bed

cut -f 1-6  loops_TNBC_unique_10000.bedpe  | tail -n+2 > loops_TNBC_unique_10000_16.bedpe

BedpeAnnotateFromBed -bed /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_SE_sort_with_SE.bed -in /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000.bedpe -out loops_TNBC_unique_10000_SEed.bedpe -col 6 -col_name SE_anno



#����TNBC�����Ժϲ�������ֻ����unique��loop�ļ�
cut -f 1-6  loops_TNBC_unique_10000.bedpe  | tail -n+2 > loops_TNBC_unique_10000_16.bedpe
bedtools intersect -a /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000_16.bedpe  -b /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_SE_sort.bed  /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_TE_sort.bed  /mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19_promoter_sort.bed -wa -wb -loj -sortout -names SE TE P
                       
                       
                       
                       

��������loop�ļ�������һ����ʶ��
awk 'BEGIN{OFS="\t"} {print $0, "loop"NR-1}' /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000.bedpe > loops_TNBC_unique_10000_withid.bedpe
����ʹ����ʱ�ļ���ת�����ôδ�ת��:���ٶ���ÿһ��bed�ļ�����ֻȡ��Ӧ��bed3��ʽ�У�Ȼ��һ���һ�ж�ȥ��
awk 'BEGIN{OFS="\t"} {print $0, "loop"NR-1}' /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000.bedpe >  /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000_withid.bedpe
cut -f 1-3,12 /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000_withid.bedpe | tail -n+2 > temp.bed
bedtools intersect -a temp.bed -b /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_SE_sort.bed /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_TE_sort.bed /mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19_promoter_sort.bed -wa -wb -loj -sortout -names SE TE P
rm temp.bed

cut -f 4-6,12 /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/TNBC_preprocessing_any/loops_TNBC_unique_10000_withid.bedpe | tail -n+2 > temp.bed
bedtools intersect -a temp.bed -b /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_SE_sort.bed /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed/TNBC_TE_sort.bed /mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19_promoter_sort.bed -wa -wb -loj -sortout -names SE TE P
rm temp.bed

������һ�����⣬������ν����ߵ����ݺ���һ�𣬱��밴��loop��Ե�˳��

#������3�߷ֿ���������ȡ����SE���ص�gene
Ҳ���Զ�HMEC��һ�����������TNBC��SE��EPC-HMEC��unique�е�SE-PC��gene��Ϊ���ģ������



