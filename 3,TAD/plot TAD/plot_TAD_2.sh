#!/bin/bash

Region_for_comparison=$1
Folder=$2
file1=$3
file2=$4


#mamba activate HiCExplorer
hicPlotTADs --tracks TADs-plotsettings.ini --region ${Region_for_comparison}  -o TAD_calling_${Region_for_comparison}_plot.png

#chr18:3000000-5500000,例如此类格式,此处区域不定，可以按照之前绘制的RCP-decay图等来指定染色体，比如我看decay图上chr22后期Mb差异比较大，就选了     
#--------------------------


hicPlotMatrix --matrix ${Folder}/${file1}  --region ${Region_for_comparison} --title ${file1}_${Region_for_comparison} --clearMaskedBins --dpi 300 --outFileName ${file1}_${Region_for_comparison}_plot.png

hicPlotMatrix --matrix ${Folder}/${file2} --region ${Region_for_comparison} --title ${file2}_${Region_for_comparison} --clearMaskedBins --dpi 300 --outFileName ${file2}_${Region_for_comparison}_plot.png

#因为每次报错都是针对我的第2个矩阵文件，所以应该是只能输入一个文件，所以分开来绘制在上面plot中的比对区域以及文件
#目前只能按照file1,2来指定对应的矩阵绘制了,都是原始的h5文件，注意和ini文件中对应