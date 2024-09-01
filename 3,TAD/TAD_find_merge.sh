#!/bin/bash

#mamba activate HiCExplorer
export HDF5_USE_FILE_LOCKING='FALSE'

#Resolution=$1           #先使用10,20,40kb的数据

#--------------------------
#原先使用h5文件    
#for file in *.h5; do
  # prefix = $(basename ${file} .h5) 
#  hicFindTADs -m ${file} --outPrefix  $(basename ${file} .h5) --minDepth $((3*${Resolution})) --maxDepth $((10*${Resolution})) --step ${Resolution} --thresholdComparisons 0.01 --delta 0.01 --correctForMultipleTesting fdr -p 15


#folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/norm_correct_mcool
#myres=(40000 20000)  #call TAD还是只用2个分辨率即可
#myres=(10000 100000 500000)   #merge domain需要，另外call 10kb+100kb+500kb的数据，最后可以merge在一起
#myrep=(HMEC_rep1 HMEC_rep2 BT549_rep1 BT549_rep2 HCC70_rep1 HCC70_rep2 MB231)
#for res in ${myres[@]};do
#  for rep in ${myrep[@]};do
#    hicFindTADs -m ${folder}/${rep}.mcool::/resolutions/${res} --outPrefix ${res}_${rep} --minDepth $((3*${res})) --maxDepth $((10*${res})) --step ${res} --thresholdComparisons 0.01 --delta 0.01 --correctForMultipleTesting fdr -p 20
#    done;
#done
    
#step还是修改为res,因为之前报错起码得是bin size


#1，call TAD,现在使用hic-pro上游合并fq之后的数据，还是使用10+20+40+100kb左右的数据
folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/merge_fq_now
myres=(10000 20000 40000 100000)
myrep=(HMEC BT549 HCC70 MB231)
for rep in ${myrep[@]};do
  for res in ${myres[@]};do
    hicFindTADs -m ${folder}/${rep}.mcool::/resolutions/${res} --outPrefix ${res}_${rep} --minDepth $((3*${res})) --maxDepth $((10*${res})) --step ${res} --thresholdComparisons 0.01 --delta 0.01 --correctForMultipleTesting fdr -p 20
  done;
  #2,然后就是merge TAD操作，注意是在同一个rep内部，合并10-20-40-100kb的res数据，注意上面res操作完成之后会获得res_rep_domains.bed文件
  hicMergeDomains --domainFiles 10000_${rep}_domains.bed 20000_${rep}_domains.bed 40000_${rep}_domains.bed 100000_${rep}_domains.bed   --outputMergedList  ${rep}_merged_noCTCF --outputRelationList  ${rep}_relation_noCTCF  --outputTreePlotPrefix ${rep}_plot_noCTCF  --outputTreePlotFormat pdf  
done

#注意上面并没有获取对应的CTCF的数据，可以之后合并的时候再考虑一下顺便使用一下500kb左右的数据

    
