# 读取左右anchor的注释数据
left_anchor_data = {}
right_anchor_data = {}

with open("MB231_MB231_leftanchor.bed", "r") as left_file:
    for line in left_file:
        parts = line.strip().split("\t")
        loop_id = parts[3]
        if loop_id not in left_anchor_data:
            left_anchor_data[loop_id] = []
        left_anchor_data[loop_id].append(parts)

with open("MB231_MB231_rightanchor.bed", "r") as right_file:
    for line in right_file:
        parts = line.strip().split("\t")
        loop_id = parts[3]
        if loop_id not in right_anchor_data:
            right_anchor_data[loop_id] = []
        right_anchor_data[loop_id].append(parts)

# 将左右anchor相同loop的数据组合放在同一行
combined_data = []

for loop_id, left_data in left_anchor_data.items():
    if loop_id in right_anchor_data:
        right_data = right_anchor_data[loop_id]
        for left_row in left_data:
            for right_row in right_data:
                combined_row = left_row + right_row[4:]  # 合并两行数据
                combined_data.append(combined_row)

# 输出组合后的数据到文件
with open("MB231_MB231_loop_annotation.bed", "w") as output_file:
    for row in combined_data:
        output_file.write("\t".join(row) + "\n")
        
        
#然后下面是对HMEC相应
left_anchor_data = {}
right_anchor_data = {}

with open("MB231_HMEC_leftanchor.bed", "r") as left_file:
    for line in left_file:
        parts = line.strip().split("\t")
        loop_id = parts[3]
        if loop_id not in left_anchor_data:
            left_anchor_data[loop_id] = []
        left_anchor_data[loop_id].append(parts)

with open("MB231_HMEC_rightanchor.bed", "r") as right_file:
    for line in right_file:
        parts = line.strip().split("\t")
        loop_id = parts[3]
        if loop_id not in right_anchor_data:
            right_anchor_data[loop_id] = []
        right_anchor_data[loop_id].append(parts)

# 将左右anchor相同loop的数据组合放在同一行
combined_data = []

for loop_id, left_data in left_anchor_data.items():
    if loop_id in right_anchor_data:
        right_data = right_anchor_data[loop_id]
        for left_row in left_data:
            for right_row in right_data:
                combined_row = left_row + right_row[4:]  # 合并两行数据
                combined_data.append(combined_row)

# 输出组合后的数据到文件
with open("MB231_HMEC_loop_annotation.bed", "w") as output_file:
    for row in combined_data:
        output_file.write("\t".join(row) + "\n")
