# -*- coding: utf-8 -*-

# 读取文件内容并拼接
file_names = ["BT549_10000.tsv", "HCC70_10000.tsv", "MB231_10000.tsv"]
output_file = "merged_file.tsv"

with open(output_file, "w") as outfile:
    for filename in file_names:
        with open(filename, "r") as infile:
            # 跳过每个文件的第一行，即列名行
            next(infile)
            for line in infile:
                outfile.write(line)

print("The files were concatenated and the results were saved in", output_file)
