#!/bin/bash

#暂时就以HMEC以及TNBC中各自的SE作为比对，看看效果
loopdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result
allEdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed
Efile=(SE TE allE)
#myrep=(TNBC HMEC) #HMEC BT549 HCC70 MB231 TNBC
declare -A suffix=(
    ["HMEC"]="HMEC_unique HMEC_common HMEC_all" 
    ["TNBC"]="TNBC_unique TNBC_common TNBC_all"
)


#for TNBC,还是分文件夹执行操作
myrep=(TNBC HMEC) #BT549 HCC70 MB231 TNBC同理
mkdir -p ${loopdir}/TNBC_preprocessing_any/EPC
for rep in ${myrep[@]};do
  for E in ${Efile[@]};do
    fs_array=(${suffix[$rep]})
    for fs in  ${fs_array[@]};do
      python /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/EPC/find_functional_loops.py ${allEdir}/${rep}_${E}_sort.bed  ${loopdir}/TNBC_preprocessing_any/loops_${fs}_10000_summit.bed  ${loopdir}/TNBC_preprocessing_any/EPC/${rep}_${E}PC_${fs}.txt
    done;
  done;
done

#后面可以更改上面的suffix后缀，HMEC不动，但是3个TNBC都设置，然后rep的话就设置各自的比对中，然后注意loop文件以及E文件地址的设置，最后就是输出，可以放在对应的文件夹下面的EPC