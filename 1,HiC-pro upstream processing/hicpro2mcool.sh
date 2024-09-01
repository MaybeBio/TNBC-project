#!/bin/bash

#  mamba activate hicexplorer
#   export HDF5_USE_FILE_LOCKING='FALSE'

#--------------------------
FOLDER=$1   #指定/mnt/disk4/haitao/bysj_seu/geo_data/hic/hicpro_output/
Resolution=$2        #指定分辨率，按照原来脚本，还是先处理10kb的
 
for folder in $FOLDER/logs/*/; do

	SAMPLE=$(basename $folder)
  hicproFILE_=/mnt/disk4/haitao/bysj_seu/geo_data/hic/hicpro_output/hic_results/matrix/${SAMPLE}/raw/${Resolution}/${SAMPLE}_${Resolution}.matrix    #暂时还是使用raw的格式，配套
  hicproBEDFILE_=/mnt/disk4/haitao/bysj_seu/geo_data/hic/hicpro_output/hic_results/matrix/${SAMPLE}/raw/${Resolution}/${SAMPLE}_${Resolution}_abs.bed

#--------------------------

  fname=$(basename ${hicproFILE_} .matrix)

  convertedFILE_=$fname".cool"

#--------------------------
# Convert HiC-Pro matrices to h5 matrices
#--------------------------

  hicConvertFormat --matrices ${hicproFILE_} --bedFileHicpro ${hicproBEDFILE_} --outFileName ${convertedFILE_} --inputFormat hicpro --outputFormat cool 

done