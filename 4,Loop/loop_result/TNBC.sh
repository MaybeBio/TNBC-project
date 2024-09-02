#!/bin/bash
#tail -n +1 BT549_10000.tsv > combined.tsv
#tail -n +2 HCC70_10000.tsv >> combined.tsv
#tail -n +2 MB231_10000.tsv >> combined.tsv

awk '{
    key = $1FS$2FS$3FS$4FS$5FS$6; # 使用前6列作为键值
    if (key in seen) {
        for (i = 7; i <= NF; i++) {
            sum[i] += $i; # 求和
        }
        count[key]++; # 计数
        next; # 跳过重复行
    }
    seen[key] = 1; # 标记为已见过的行
    print $0; # 输出第一次出现的行
}
END {
    for (key in count) {
        for (i = 7; i <= NF; i++) {
            printf("%s%s", sum[i]/count[key], (i<NF) ? FS : RS); # 计算均值并输出
        }
    }
}' combined.tsv > combined_processed.tsv
