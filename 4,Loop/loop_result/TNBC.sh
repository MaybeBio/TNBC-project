#!/bin/bash
#tail -n +1 BT549_10000.tsv > combined.tsv
#tail -n +2 HCC70_10000.tsv >> combined.tsv
#tail -n +2 MB231_10000.tsv >> combined.tsv

awk '{
    key = $1FS$2FS$3FS$4FS$5FS$6; # ʹ��ǰ6����Ϊ��ֵ
    if (key in seen) {
        for (i = 7; i <= NF; i++) {
            sum[i] += $i; # ���
        }
        count[key]++; # ����
        next; # �����ظ���
    }
    seen[key] = 1; # ���Ϊ�Ѽ�������
    print $0; # �����һ�γ��ֵ���
}
END {
    for (key in count) {
        for (i = 7; i <= NF; i++) {
            printf("%s%s", sum[i]/count[key], (i<NF) ? FS : RS); # �����ֵ�����
        }
    }
}' combined.tsv > combined_processed.tsv
