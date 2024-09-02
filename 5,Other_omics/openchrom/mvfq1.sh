#!/bin/bash

# 创建fq1文件夹（如果不存在）
mkdir -p fq1

# 循环处理匹配的文件
for file in *_val*.fq; do
    # 获取文件名和后缀
    filename=$(basename "$file")
    extension="${filename##*.}"

    # 检查文件名中_val的出现次数
    val_count=$(echo "$filename" | grep -o '_val' | wc -l)

    # 如果_val出现一次，则移动文件到fq1文件夹中
    if [ "$val_count" -eq 1 ]; then
        mv "$file" "fq1/$filename"
        echo "已移动文件: $file 到 fq1 文件夹"
    fi
done

