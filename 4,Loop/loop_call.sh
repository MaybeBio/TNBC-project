#!/bin/bash

# mamba activate mustache

folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/merge_fq_now
myrep=(HMEC BT549 HCC70 MB231) 
myres=(10000 20000 40000) #��ʱ��ʹ����3��res�����ݣ�������ݽű��ϲ���һ�����ֲ�res�Ĳ��㣬��ʵ�ǽ���ʹ��10+20kbΪ��

for rep in ${myrep[@]};do
  for res in ${myres[@]};do
    python3 /mnt/disk4/haitao/software/mustache/mustache/mustache.py -f ${folder}/${rep}.mcool -r ${res} -pt 0.05 -p 15 -o ${rep}_${res}.tsv 
  done;
done
#�˴�������ʱ�����norm����Ϊ�ᱨ���ο�https://github.com/ay-lab/mustache/issues/26
#�������res���ݣ�����ʹ��������ʽ����ָ��