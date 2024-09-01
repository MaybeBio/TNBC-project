#!/bin/bash


echo ${SAMPLE}

#mamba activate HiCExplorer
export HDF5_USE_FILE_LOCKING='FALSE'

mkdir -p ${SAMPLE}
cd ${SAMPLE}

Resolution=$1           #先使用10,20,40kb的数据

#--------------------------
cd /mnt/disk4/haitao/bysj_seu/geo_data/hic/script2/${Resolution}

for file in *.h5; do
        prefix=$(basename "$file" .h5)

#--------------------------
# TADs detection
#--------------------------

        hicFindTADs -m ${file} \       
        --outPrefix ${prefix} \
        --minDepth 3*${Resolution} \
        --maxDepth 10*${Resolution} \
        --step 10000 \
        --thresholdComparisons 0.01 \
        --delta 0.01 \
        --correctForMultipleTesting fdr \
        -p 20
#矩阵用的还是矫正之后的，可以用iced的
#输出文件可以使用hicPlotTADs绘画
#--minDepth至少是res的3倍，所以我如果用10kb来识别TAD，就可以使用30kb；但是我如果使用20kb来识别的话，就是60kb；
#--maxDepth至少是res的6-10倍，所以同上是60kb-100kb，或者120kb-2000kb,目前就使用10、20，40
#--thresholdComparisons脚本中是0.05，此处使用默认的0.01

#下面这一段其实是可以直接分离的

#--------------------------
# Plot TADs (Matplotlib must be > 3.1.1)
#--------------------------

hicPlotTADs --tracks TADs-plotsettings.ini --region chr18:3000000-5500000 -o TAD_calling_plot.png

#--------------------------

hicPlotMatrix -m ${FILE_} \
        --region chr18:1000000-8000000 \
        --title ${SAMPLE} \
        --clearMaskedBins \
        --dpi 300 \
        --outFileName chr18_plot.png
