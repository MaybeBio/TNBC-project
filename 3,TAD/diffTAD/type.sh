#!/bin/bash

#mamba activate HiCExplorer

diffTADfolder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script4/diffTAD
TADfolder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script4/merge_fq
myrep=(BT549 HCC70 MB231) # HMEC作为control
myres=(40000 20000 10000)

for res in ${myres[@]};do
  for rep in ${myrep[@]};do
  mkdir -p ${res}_${rep}_diff
    python /mnt/disk4/haitao/software/diff_tad_type/diff_tad_type.py --diff_tad_file ${diffTADfolder}/${res}_HMEC_vs_${rep}_rejected.diff_tad  --test_tad_file ${TADfolder}/${res}_${rep}_domains.bed  --control_tad_file ${TADfolder}/${res}_HMEC_domains.bed  --out_dir ${res}_${rep}_diff 
   done;
done  
    
#percent 0.1就是原文中的percent