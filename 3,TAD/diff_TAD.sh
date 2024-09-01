#!/bin/bash

#mamba activate HiCExplorer

#使用hicDifferentialTAD，因为使用的mcool文件已经经过了norm+correct，可以正常使用
myrep=(BT549 HCC70 MB231) # HMEC作为control
myres=(40000 20000 10000) #暂时先使用20kb以及40kb的数据
folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/merge_fq_now
matrixfolder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script4/merge_fq

for res in ${myres[@]};do
  for rep in ${myrep[@]};do
    hicDifferentialTAD -tm  ${folder}/${rep}.mcool::/resolutions/${res} -cm ${folder}/HMEC.mcool::/resolutions/${res} -td  ${matrixfolder}/${res}_${rep}_domains.bed   -o ${res}_HMEC_vs_${rep}   -p 0.05 -m all -mr all -t 15
  done;
done
    

#control vs treatment,需要的是treatment的TAD文件+hic的cool文件，control的hic的cool文件;另外mr还是和m一致，全all或全one