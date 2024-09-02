#!/bin/bash

# 遍历当前目录下所有的txt文件
# 从.txt文件生成.bed文件并按需分割字段
for file in *.txt; do
    prefix=$(echo "$file" | cut -d'_' -f1)

    if [ "$prefix" == "HCC70" ] || [ "$prefix" == "BT549" ]; then
        awk 'NR>5 {print $2 "\t" $3 "\t" $4 "\t" $6 "\t" $9}' FS='\t' OFS='\t' "$file" > "${file%.txt}.bed"
    else
        awk 'NR>5 {print $2 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' FS='\t' OFS='\t' "$file" > "${file%.txt}.bed"
    fi
done


mv BT549_peaks_AllStitched.table.bed BT549_allE.bed
mv BT549_peaks_SuperStitched.table.bed BT549_SE.bed
mv HCC70_peaks_AllStitched.table.bed HCC70_allE.bed
mv HCC70_peaks_SuperStitched.table.bed HCC70_SE.bed
mv HMEC_H3K27ac_peaks_AllStitched.table.bed  HMEC_allE.bed
mv HMEC_H3K27ac_peaks_SuperStitched.table.bed HMEC_SE.bed
mv MB231_H3K27ac_peaks_AllStitched.table.bed MB231_allE.bed
mv MB231_H3K27ac_peaks_SuperStitched.table.bed MB231_SE.bed