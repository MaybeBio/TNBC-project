#!/bin/bash

#mamba activate hicexplorer
Folder=$1 #ָ��/mnt/disk4/haitao/bysj_seu/geo_data/hic/hic2juice/
#�����ʹ��geo��ԭʼhic����ʵ��Ļ���ָ���ļ���~/bysj_seu/geo_data/hic/raw_geo_data

for file in ${Folder}/*.hic; do

    prefix=$(basename ${file} .hic)  #������ֵ�����пո񣬿��Կ�������ɫ��֮ǰ�Ľű�Ҳһ��
    #hicConvertFormat -m ${file} --inputFormat hic --outputFormat cool -o ${prefix}.mcool  --resolutions  20000 40000
    
    hic2cool convert ${file} ${prefix}_.mcool -r 0 -p 20
done






