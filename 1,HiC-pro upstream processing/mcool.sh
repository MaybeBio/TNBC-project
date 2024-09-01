#!/bin/bash

#mamba activate hicexplorer
Folder=$1 #指定/mnt/disk4/haitao/bysj_seu/geo_data/hic/hic2juice/
#如果是使用geo中原始hic数据实验的话就指定文件夹~/bysj_seu/geo_data/hic/raw_geo_data

for file in ${Folder}/*.hic; do

    prefix=$(basename ${file} .hic)  #变量赋值不能有空格，可以看字体颜色，之前的脚本也一样
    #hicConvertFormat -m ${file} --inputFormat hic --outputFormat cool -o ${prefix}.mcool  --resolutions  20000 40000
    
    hic2cool convert ${file} ${prefix}_.mcool -r 0 -p 20
done






