#!/bin/bash

# mamba activate mustache

folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/merge_fq_now
myrep=(HMEC BT549 HCC70 MB231) 
myres=(10000 20000 40000) #暂时就使用这3个res的数据，最后依据脚本合并在一起，以弥补res的不足，其实是建议使用10+20kb为主

for rep in ${myrep[@]};do
  for res in ${myres[@]};do
    python3 /mnt/disk4/haitao/software/mustache/mustache/mustache.py -f ${folder}/${rep}.mcool -r ${res} -pt 0.05 -p 15 -o ${rep}_${res}.tsv 
  done;
done
#此处还是暂时不添加norm，因为会报错，参考https://github.com/ay-lab/mustache/issues/26
#另外对于res数据，还是使用数字形式进行指定