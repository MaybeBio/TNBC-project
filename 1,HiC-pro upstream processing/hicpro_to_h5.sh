#!/bin/bash

#  mamba activate hicexplorer
#   export HDF5_USE_FILE_LOCKING='FALSE'

#--------------------------
FOLDER=$1   #ָ��/mnt/disk4/haitao/bysj_seu/geo_data/hic/hicpro_output/
Resolution=$2        #ָ���ֱ��ʣ�����ԭ���ű��������ȴ���10kb��
 
for folder in $FOLDER/logs/*/; do

	SAMPLE=$(basename $folder)
  hicproFILE_=/mnt/disk4/haitao/bysj_seu/geo_data/hic/hicpro_output/hic_results/matrix/${SAMPLE}/iced/${Resolution}/${SAMPLE}_${Resolution}_iced.matrix    #��ʱ���������ָ���ֱ��ʣ���һ����10kb��������������һ��ѭ������ο���15��
  hicproBEDFILE_=/mnt/disk4/haitao/bysj_seu/geo_data/hic/hicpro_output/hic_results/matrix/${SAMPLE}/raw/${Resolution}/${SAMPLE}_${Resolution}_abs.bed

#--------------------------

  fname=$(basename ${hicproFILE_} .matrix)

  convertedFILE_=$fname".h5"

#--------------------------
# Convert HiC-Pro matrices to h5 matrices
#--------------------------

  hicConvertFormat --matrices ${hicproFILE_} --outFileName ${convertedFILE_} --inputFormat hicpro --outputFormat h5 --bedFileHicpro ${hicproBEDFILE_}

done