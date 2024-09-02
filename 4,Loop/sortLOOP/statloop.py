print("TNBC_TNBC")
# 读取注释文件并提取所需的列
annotations = []

with open("TNBC_TNBC_loop_annotation.bed", "r") as annotation_file:
    for line in annotation_file:
        parts = line.strip().split("\t")
        annotations.append((parts[4], parts[8]))

# 统计组合的频次
combination_counts = {}

for annotation_pair in annotations:
    if annotation_pair not in combination_counts:
        combination_counts[annotation_pair] = 0
    combination_counts[annotation_pair] += 1

# 输出统计结果
for combination, count in combination_counts.items():
    print(f"{combination[0]} - {combination[1]}: {count}")

print("TNBC_HMEC")
# 读取注释文件并提取所需的列
annotations = []

with open("TNBC_HMEC_loop_annotation.bed", "r") as annotation_file:
    for line in annotation_file:
        parts = line.strip().split("\t")
        annotations.append((parts[4], parts[8]))

# 统计组合的频次
combination_counts = {}

for annotation_pair in annotations:
    if annotation_pair not in combination_counts:
        combination_counts[annotation_pair] = 0
    combination_counts[annotation_pair] += 1

# 输出统计结果
for combination, count in combination_counts.items():
    print(f"{combination[0]} - {combination[1]}: {count}")