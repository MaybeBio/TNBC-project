#!/bin/bash


echo ${SAMPLE}

#mamba activate HiCExplorer
export HDF5_USE_FILE_LOCKING='FALSE'

mkdir -p ${SAMPLE}
cd ${SAMPLE}

Resolution=$1           #��ʹ��10,20,40kb������

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
#�����õĻ��ǽ���֮��ģ�������iced��
#����ļ�����ʹ��hicPlotTADs�滭
#--minDepth������res��3���������������10kb��ʶ��TAD���Ϳ���ʹ��30kb�����������ʹ��20kb��ʶ��Ļ�������60kb��
#--maxDepth������res��6-10��������ͬ����60kb-100kb������120kb-2000kb,Ŀǰ��ʹ��10��20��40
#--thresholdComparisons�ű�����0.05���˴�ʹ��Ĭ�ϵ�0.01

#������һ����ʵ�ǿ���ֱ�ӷ����

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
