#!/bin/bash

#mamba activate HiCExplorer
export HDF5_USE_FILE_LOCKING='FALSE'

#folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/norm_correct_mcool
#myres=(100000 500000)
#for res in ${myres[@]};do 
  #hicPlotDistVsCounts --plotFile HMEC_vs_TNBC_all_${res}.png -m ${folder}/HMEC_rep2.mcool::/resolutions/${res} ${folder}/HMEC_rep1.mcool::/resolutions/${res} ${folder}/HCC70_rep2.mcool::/resolutions/${res} ${folder}/HCC70_rep1.mcool::/resolutions/${res} ${folder}/BT549_rep2.mcool::/resolutions/${res} ${folder}/BT549_rep1.mcool::/resolutions/${res} ${folder}/MB231.mcool::/resolutions/${res} --labels 'HMEC_rep2' 'HMEC_rep1' 'HCC70_rep2' 'HCC70_rep1' 'BT549_rep2' 'BT549_rep1' 'MB231' --outFileData HMEC_vs_TNBC_all_${res}.txt --plotsize 6 5
  #hicPlotDistVsCounts --plotFile HMEC_vs_TNBC_perchr_${res}.png -m ${folder}/HMEC_rep2.mcool::/resolutions/${res} ${folder}/HMEC_rep1.mcool::/resolutions/${res} ${folder}/HCC70_rep2.mcool::/resolutions/${res} ${folder}/HCC70_rep1.mcool::/resolutions/${res} ${folder}/BT549_rep2.mcool::/resolutions/${res} ${folder}/BT549_rep1.mcool::/resolutions/${res} ${folder}/MB231.mcool::/resolutions/${res} --labels 'HMEC_rep2' 'HMEC_rep1' 'HCC70_rep2' 'HCC70_rep1' 'BT549_rep2' 'BT549_rep1' 'MB231' --outFileData HMEC_vs_TNBC_perchr_${res}.txt --plotsize 6 5 --perchr
#done

#上面是分rep的，下面是合并之后的
folder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/merge_fq_now
myres=(100000 500000)
for res in ${myres[@]};do 
  hicPlotDistVsCounts --plotFile HMEC_vs_TNBC_allmerge_${res}.png -m ${folder}/HMEC.mcool::/resolutions/${res}  ${folder}/HCC70.mcool::/resolutions/${res}  ${folder}/BT549.mcool::/resolutions/${res}  ${folder}/MB231.mcool::/resolutions/${res} --labels 'HMEC' 'HCC70' 'BT549' 'MB231' --outFileData HMEC_vs_TNBC_allmerge_${res}.txt --plotsize 6 5
  hicPlotDistVsCounts --plotFile HMEC_vs_TNBC_perchrmerge_${res}.png -m ${folder}/HMEC.mcool::/resolutions/${res}  ${folder}/HCC70.mcool::/resolutions/${res}  ${folder}/BT549.mcool::/resolutions/${res}  ${folder}/MB231.mcool::/resolutions/${res} --labels 'HMEC' 'HCC70' 'BT549' 'MB231' --outFileData HMEC_vs_TNBC_perchrmerge_${res}.txt --plotsize 6 5 --perchr
done