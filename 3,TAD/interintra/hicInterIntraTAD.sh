#!/bin/bash

#mamba activate HiCExplorer
myrep=(HMEC BT549 HCC70 MB231) 
myres=(10000 20000 40000) #暂时先使用20kb以及40kb的数据
folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/merge_fq_now
matrixfolder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script4/merge_fq

for res in ${myres[@]};do
  for rep in ${myrep[@]};do
    hicInterIntraTAD -m ${folder}/${rep}.mcool::/resolutions/${res} -td ${matrixfolder}/${res}_${rep}_domains.bed -o ${res}_${rep} -op ${res}_${rep} -t 15 
  done;
done