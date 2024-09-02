#!/bin/bash

mcoolfolder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/merge_fq_now
myrep=(BT549 HCC70 MB231)

for rep in ${myrep[@]};do
  mkdir -p ${rep}
  python3 /mnt/disk4/haitao/software/mustache/mustache/diff_mustache.py -f1 ${mcoolfolder}/${rep}.mcool -f2 ${mcoolfolder}/HMEC.mcool  -pt 0.05 -pt2 0.05  -p 15 -o ${rep} -r 10000         
done