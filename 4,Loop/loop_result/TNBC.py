# -*- coding: utf-8 -*-

# ��ȡ�ļ����ݲ�ƴ��
file_names = ["BT549_10000.tsv", "HCC70_10000.tsv", "MB231_10000.tsv"]
output_file = "merged_file.tsv"

with open(output_file, "w") as outfile:
    for filename in file_names:
        with open(filename, "r") as infile:
            # ����ÿ���ļ��ĵ�һ�У���������
            next(infile)
            for line in infile:
                outfile.write(line)

print("The files were concatenated and the results were saved in", output_file)
