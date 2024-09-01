#!/bin/bash

#mamba activate HiCExplorer

region=$2    #指定绘制的区域，类似于chrX:6800000-8500000，格式和bed文件中类似以及和一些风险位点也类似
ini=$1   #指定绘制的ini配置文件
output_prefix=$3  #指定输出图片的前缀

#common比较，即全局TAD track图景比较，不涉及特异性TAD，仅仅是show，使用TADs-common.ini，至于差异的TAD比较可以另外设置一个ini文件
#可以在ini中将所有的TAD都统一起来，也可以写个循环，然后设置不同的ini，到时候在循环里统一绘制同一个区域

hicPlotTADs --tracks  ${ini} --region ${region} -o ${output_prefix}.pdf