#!/bin/bash
# Script to use https://github.com/magmarshh/TADWidth?tab=readme-ov-file 
#使用hicexplorer中hicfindTAD识别的TAD的domain的bed文件，使用40kb+20kb的数据，10kb的就不用了

TADdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script4/merge_fq/

myres=(40000 20000)
myrep=(HMEC BT549 HCC70 MB231)
for res in ${myres[@]};do
  for rep in ${myrep[@]};do
    python /mnt/disk4/haitao/software/TADWidth/TADWidth_piechart.py --tad ${TADdir}/${res}_${rep}_domains.bed  --bins 10Kb-50Kb+ --output ${res}_${rep}_sizes_piechart.pdf
    python /mnt/disk4/haitao/software/TADWidth/TADWidth_piechart.py --tad ${TADdir}/10000_${rep}_domains.bed  --res 10000 --output 10000_${rep}_sizes_piechart.pdf
    python /mnt/disk4/haitao/software/TADWidth/TADWidth_violinplot.py --tad ${TADdir}/${res}_${rep}_domains.bed --labels ${rep} --output ${res}_${rep}_sizes_violinplot.pdf
  done
    python /mnt/disk4/haitao/software/TADWidth/TADWidth_violinplot.py --tad ${TADdir}/${res}_HMEC_domains.bed,${TADdir}/${res}_BT549_domains.bed,${TADdir}/${res}_HCC70_domains.bed,${TADdir}/${res}_MB231_domains.bed --labels HMEC,BT549,HCC70,MB231    --figWidth 10 --figHeight 7 --output ${res}_sizes_violinplot.pdf

done



