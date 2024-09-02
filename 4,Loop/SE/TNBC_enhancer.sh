#!/bin/bash

#对所有的allE文件进行提取，提取
for file in *_allE.bed; do
    head -n 1 "$file" > "${file%_allE.bed}_TE.bed"
    awk '$5 == 0 {print}' "$file" >> "${file%_allE.bed}_TE.bed"
done


#TE
tail -n +1 "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/BT549_TE.bed" > TNBC_TE.bed
tail -n +2 "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/HCC70_TE.bed" >> TNBC_TE.bed
tail -n +2 "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/MB231_TE.bed" >> TNBC_TE.bed

#SE
tail -n +1 "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/BT549_SE.bed" > TNBC_SE.bed
tail -n +2 "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/HCC70_SE.bed" >> TNBC_SE.bed
tail -n +2 "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/MB231_SE.bed" >> TNBC_SE.bed

#ALL E
tail -n +1 "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/BT549_allE.bed" > TNBC_allE.bed
tail -n +2 "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/HCC70_allE.bed" >> TNBC_allE.bed
tail -n +2 "/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/MB231_allE.bed" >> TNBC_allE.bed