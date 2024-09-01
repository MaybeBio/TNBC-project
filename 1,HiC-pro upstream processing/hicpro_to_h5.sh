#!/bin/bash

#  mamba activate hicexplorer
#   export HDF5_USE_FILE_LOCKING='FALSE'

#--------------------------
FOLDER=$1   #指定/mnt/disk4/haitao/bysj_seu/geo_data/hic/hicpro_output/
Resolution=$2        #指定分辨率，按照原来脚本，还是先处理10kb的
 
for folder in $FOLDER/logs/*/; do

	SAMPLE=$(basename $folder)
  hicproFILE_=/mnt/disk4/haitao/bysj_seu/geo_data/hic/hicpro_output/hic_results/matrix/${SAMPLE}/iced/${Resolution}/${SAMPLE}_${Resolution}_iced.matrix    #到时候这里可以指定分辨率，不一定得10kb；或者外面再套一层循环；或参考第15行
  hicproBEDFILE_=/mnt/disk4/haitao/bysj_seu/geo_data/hic/hicpro_output/hic_results/matrix/${SAMPLE}/raw/${Resolution}/${SAMPLE}_${Resolution}_abs.bed

#--------------------------

  fname=$(basename ${hicproFILE_} .matrix)

  convertedFILE_=$fname".h5"

#--------------------------
# Convert HiC-Pro matrices to h5 matrices
#--------------------------

  hicConvertFormat --matrices ${hicproFILE_} --outFileName ${convertedFILE_} --inputFormat hicpro --outputFormat h5 --bedFileHicpro ${hicproBEDFILE_}

done