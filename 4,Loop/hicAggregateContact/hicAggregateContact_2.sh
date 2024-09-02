#!/bin/bash

#mamba activate hicexplorer

myrep=(HMEC BT549 HCC70 MB231)
mcooldir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/merge_fq_now #mcool file
myres=(40000 20000)       #按GENOVA中的PESCAn的res的需求，建议使用的res是20kb-40kb左右,PE-SCAn 20kb-40kb
plottype=(2d 3d)
#以下输入的区域基本上都是以HMEC为主 
bed1=$1 #输入比较的第一个区域的bed文件
bed2=$2 #输入比较的第二个区域的bed文件
 
for res in ${myres[@]};do
  for rep in ${myrep[@]};do
    for plot in ${plottype[@]};do
  
  base1=$(basename ${bed1})
  base2=$(basename ${bed2})
  prefix_bed1=${base1%.*}
  prefix_bed2=${base2%.*}
  
    hicAggregateContacts  -m ${mcooldir}/${rep}.mcool::/resolutions/${res} --outFileName ${res}_${rep}_${prefix_bed1}_${prefix_bed2}_${plot} --BED ${bed1} --BED2 ${bed2} --transform obs/exp --operationType mean --mode intra-chr --largeRegionsOperation center --plotType ${plot}  
    
    done;
  done;
done
    
#--considerStrandDirection 这个参数在CTCF的时候考虑使用
#--outFilePrefixMatrix  --outFileContactPairs --outFileObsExp --diagnosticHeatmapFile 后面的这些选项就全都暂时不使用，后续有需求再深入
#总之使用两个bed文件的可以和GENOVA中的CSCAn工具互补，然后使用一个bed文件的可以和GENOVA中的PESCAn工具互补
