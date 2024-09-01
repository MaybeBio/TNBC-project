#!/bin/bash

Region_for_comparison=$1
Folder=$2
file1=$3
file2=$4


#mamba activate HiCExplorer
hicPlotTADs --tracks TADs-plotsettings.ini --region ${Region_for_comparison}  -o TAD_calling_${Region_for_comparison}_plot.png

#chr18:3000000-5500000,��������ʽ,�˴����򲻶������԰���֮ǰ���Ƶ�RCP-decayͼ����ָ��Ⱦɫ�壬�����ҿ�decayͼ��chr22����Mb����Ƚϴ󣬾�ѡ��     
#--------------------------


hicPlotMatrix --matrix ${Folder}/${file1}  --region ${Region_for_comparison} --title ${file1}_${Region_for_comparison} --clearMaskedBins --dpi 300 --outFileName ${file1}_${Region_for_comparison}_plot.png

hicPlotMatrix --matrix ${Folder}/${file2} --region ${Region_for_comparison} --title ${file2}_${Region_for_comparison} --clearMaskedBins --dpi 300 --outFileName ${file2}_${Region_for_comparison}_plot.png

#��Ϊÿ�α���������ҵĵ�2�������ļ�������Ӧ����ֻ������һ���ļ������Էֿ�������������plot�еıȶ������Լ��ļ�
#Ŀǰֻ�ܰ���file1,2��ָ����Ӧ�ľ��������,����ԭʼ��h5�ļ���ע���ini�ļ��ж�Ӧ