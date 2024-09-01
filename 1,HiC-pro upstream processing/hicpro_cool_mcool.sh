#!/bin/bash

# mamba activate hicexplorer
# export HDF5_USE_FILE_LOCKING='FALSE'

#之前的脚本太混乱了，现在将所有脚本都整合一遍，重新跑一下，并且是分开处理

#1，先将hicpro文件转换为cool文件，都是分rep,注意以rep为主
data=/mnt/disk4/haitao/merge_hic/data_output/hic_results/matrix
myrep=(HMEC HCC70 BT549 MB231)
myres=(10000 20000 40000 100000 500000)

#下面全部使用raw的数据，还是先res处理，每一个res之后需要在同一res的基础上去norm
for res in ${myres[@]};do
  for rep in ${myrep[@]};do
    hicproFILE_=${data}/${rep}/raw/${res}/${rep}_${res}.matrix
    hicproBEDFILE_=${data}/${rep}/raw/${res}/${rep}_${res}_abs.bed
    fname=${rep}_${res}
    convertedFILE_=${fname}.cool
    #开始转化文件，获得${rep}_${res}.cool文件
    hicConvertFormat --matrices ${hicproFILE_} --bedFileHicpro ${hicproBEDFILE_} --outFileName ${convertedFILE_} --inputFormat hicpro --outputFormat cool 
    #将${rep1}_${res}.cool文件进行矫正
  done;
  
  #2，注意完成了一个res内的所有的cool文件的生成，再接着用所有文件的cool在同一res内做norm
  hicNormalize -m HMEC_${res}.cool  HCC70_${res}.cool  BT549_${res}.cool  MB231_${res}.cool --normalize smallest -o HMEC_${res}_norm.cool  HCC70_${res}_norm.cool  BT549_${res}_norm.cool  MB231_${res}_norm.cool 

  #3，注意在norm之后需要correct每一个文件，这个时候还在res循环里面，但是rep循环已经跳出，所以需要再执行rep循环，得到了${rep}_${res}_norm_KR.cool
  for rep in ${myrep[@]};do
  hicCorrectMatrix correct --matrix ${rep}_${res}_norm.cool  --correctionMethod KR --outFileName ${rep}_${res}_norm_KR.cool
  done;
done

#4，上面完成了每一个res下的rep的矫正，但是需要将最后矫正的cool文件转换为mcool文件，但是这里实际上是多res的，所以需要跳出上面的循环，重新对于每一个rep而言
for rep in ${myrep[@]};do
  hicConvertFormat -m ${rep}_10000_norm_KR.cool  ${rep}_20000_norm_KR.cool  ${rep}_40000_norm_KR.cool  ${rep}_100000_norm_KR.cool  ${rep}_500000_norm_KR.cool  --inputFormat cool --outputFormat mcool -o  ${rep}.mcool
done



