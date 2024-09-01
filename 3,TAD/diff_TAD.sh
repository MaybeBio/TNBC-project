#!/bin/bash

#mamba activate HiCExplorer

#ʹ��hicDifferentialTAD����Ϊʹ�õ�mcool�ļ��Ѿ�������norm+correct����������ʹ��
myrep=(BT549 HCC70 MB231) # HMEC��Ϊcontrol
myres=(40000 20000 10000) #��ʱ��ʹ��20kb�Լ�40kb������
folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/merge_fq_now
matrixfolder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script4/merge_fq

for res in ${myres[@]};do
  for rep in ${myrep[@]};do
    hicDifferentialTAD -tm  ${folder}/${rep}.mcool::/resolutions/${res} -cm ${folder}/HMEC.mcool::/resolutions/${res} -td  ${matrixfolder}/${res}_${rep}_domains.bed   -o ${res}_HMEC_vs_${rep}   -p 0.05 -m all -mr all -t 15
  done;
done
    

#control vs treatment,��Ҫ����treatment��TAD�ļ�+hic��cool�ļ���control��hic��cool�ļ�;����mr���Ǻ�mһ�£�ȫall��ȫone