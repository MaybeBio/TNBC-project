# 读取左右anchor的注释数据
left_anchor_data = {}
right_anchor_data = {}

with open("TNBC_TNBC_leftanchor.bed", "r") as left_file:
    for line in left_file:
        parts = line.strip().split("\t")
        loop_id = parts[3]
        if loop_id not in left_anchor_data:
            left_anchor_data[loop_id] = []
        left_anchor_data[loop_id].append(parts)

with open("TNBC_TNBC_rightanchor.bed", "r") as right_file:
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
                combined_row = left_row + right_row[4:]  # 合并两行数据，这里操作的时候左anchor的bed数据保留了下来，右anchor的bed数据没有保留，直接将feature注释合并，因为已经有了loop id作为参考，所以实际上左边anchor的bed也不需要，如果要保留完整的anchor的左右数据，可以在这里修改，将所有数据都合并上一行
                combined_data.append(combined_row)

# 输出组合后的数据到文件
with open("TNBC_TNBC_loop_annotation.bed", "w") as output_file:
    for row in combined_data:
        output_file.write("\t".join(row) + "\n")
        
        
#然后下面是对HMEC相应
left_anchor_data = {}
right_anchor_data = {}

with open("TNBC_HMEC_leftanchor.bed", "r") as left_file:
    for line in left_file:
        parts = line.strip().split("\t")
        loop_id = parts[3]
        if loop_id not in left_anchor_data:
            left_anchor_data[loop_id] = []
        left_anchor_data[loop_id].append(parts)

with open("TNBC_HMEC_rightanchor.bed", "r") as right_file:
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
with open("TNBC_HMEC_loop_annotation.bed", "w") as output_file:
    for row in combined_data:
        output_file.write("\t".join(row) + "\n")
