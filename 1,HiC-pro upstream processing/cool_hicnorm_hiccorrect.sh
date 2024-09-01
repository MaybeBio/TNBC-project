#!/bin/bash

#mamba activate hicexplorer

#cool file not cooler balance£¬better use hicexplorer 2 steps norm+correct

myrep=(BT549_rep1 BT549_rep2 HMEC_rep1 HMEC_rep2 HCC70_rep1 HCC70_rep2 MB231) 
myres=(10000 20000 40000 100000 500000)
for rep in  ${myrep[@]};do

  hicNormalize -m ${rep}.mcool::/resolutions/10000 ${rep}.mcool::/resolutions/20000 ${rep}.mcool::/resolutions/40000 ${rep}.mcool::/resolutions/100000 ${rep}.mcool::/resolutions/500000 --normalize smallest -o ${rep}_10000_norm.cool  ${rep}_20000_norm.cool  ${rep}_40000_norm.cool  ${rep}_100000_norm.cool  ${rep}_500000_norm.cool
  
  for res in ${myres[@]};do
  
  hicCorrectMatrix correct --matrix ${rep}_${res}_norm.cool  --correctionMethod KR --outFileName ${rep}_${res}_norm_KR.cool
  
  done

#can correct both in ICE and KR,now use KR 
done





  