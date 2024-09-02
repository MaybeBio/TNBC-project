# 读取两个基因集合文件
with open("/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/Fig7_HiC_Enrichments/TNBC_unique_SEPC_genes.bed") as file1:
    gene_set1 = set(file1.read().splitlines())

with open("/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/Fig7_HiC_Enrichments/HMEC_unique_SEPC_genes.bed") as file2:
    gene_set2 = set(file2.read().splitlines())

# 计算差集
gene_set1_minus_gene_set2 = gene_set1 - gene_set2
gene_set2_minus_gene_set1 = gene_set2 - gene_set1

# 输出结果
print("TNBC-HMEC:")
print(gene_set1_minus_gene_set2)

print("HMEC-TNBC:")
print(gene_set2_minus_gene_set1)

