#!/bin/bash

# mamba activate hicexplorer
# export HDF5_USE_FILE_LOCKING='FALSE'

#֮ǰ�Ľű�̫�����ˣ����ڽ����нű�������һ�飬������һ�£������Ƿֿ�����

#1���Ƚ�hicpro�ļ�ת��Ϊcool�ļ������Ƿ�rep,ע����repΪ��
data=/mnt/disk4/haitao/merge_hic/data_output/hic_results/matrix
myrep=(HMEC HCC70 BT549 MB231)
myres=(10000 20000 40000 100000 500000)

#����ȫ��ʹ��raw�����ݣ�������res����ÿһ��res֮����Ҫ��ͬһres�Ļ�����ȥnorm
for res in ${myres[@]};do
  for rep in ${myrep[@]};do
    hicproFILE_=${data}/${rep}/raw/${res}/${rep}_${res}.matrix
    hicproBEDFILE_=${data}/${rep}/raw/${res}/${rep}_${res}_abs.bed
    fname=${rep}_${res}
    convertedFILE_=${fname}.cool
    #��ʼת���ļ������${rep}_${res}.cool�ļ�
    hicConvertFormat --matrices ${hicproFILE_} --bedFileHicpro ${hicproBEDFILE_} --outFileName ${convertedFILE_} --inputFormat hicpro --outputFormat cool 
    #��${rep1}_${res}.cool�ļ����н���
  done;
  
  #2��ע�������һ��res�ڵ����е�cool�ļ������ɣ��ٽ����������ļ���cool��ͬһres����norm
  hicNormalize -m HMEC_${res}.cool  HCC70_${res}.cool  BT549_${res}.cool  MB231_${res}.cool --normalize smallest -o HMEC_${res}_norm.cool  HCC70_${res}_norm.cool  BT549_${res}_norm.cool  MB231_${res}_norm.cool 

  #3��ע����norm֮����Ҫcorrectÿһ���ļ������ʱ����resѭ�����棬����repѭ���Ѿ�������������Ҫ��ִ��repѭ�����õ���${rep}_${res}_norm_KR.cool
  for rep in ${myrep[@]};do
  hicCorrectMatrix correct --matrix ${rep}_${res}_norm.cool  --correctionMethod KR --outFileName ${rep}_${res}_norm_KR.cool
  done;
done

#4�����������ÿһ��res�µ�rep�Ľ�����������Ҫ����������cool�ļ�ת��Ϊmcool�ļ�����������ʵ�����Ƕ�res�ģ�������Ҫ���������ѭ�������¶���ÿһ��rep����
for rep in ${myrep[@]};do
  hicConvertFormat -m ${rep}_10000_norm_KR.cool  ${rep}_20000_norm_KR.cool  ${rep}_40000_norm_KR.cool  ${rep}_100000_norm_KR.cool  ${rep}_500000_norm_KR.cool  --inputFormat cool --outputFormat mcool -o  ${rep}.mcool
done



