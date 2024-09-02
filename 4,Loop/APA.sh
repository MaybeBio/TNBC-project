#!/bin/bash

myrep=(BT549 HCC70 MB231) #HMEC BT549 HCC70 MB231
hicdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/hic2juice/merge_fq #hic
loopdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result  #loop
# File suffix, corresponds to the loop conditions 
HMECsuffix=( HMEC_unique HMEC_common )
BT549suffix=( BT549_unique BT549_common )
HCC70suffix=( HCC70_unique HCC70_common )
MB231suffix=( MB231_unique MB231_common )

#上面方法行不通，还是另外自定义数组：
declare -A suffix=(
    ["BT549"]="BT549_unique BT549_common BT549_all HMEC_unique HMEC_common HMEC_all"
    ["HCC70"]="HCC70_unique HCC70_common HCC70_all HMEC_unique HMEC_common HMEC_all" 
    ["MB231"]="MB231_unique MB231_common MB231_all HMEC_unique HMEC_common HMEC_all"
)
  
for rep in ${myrep[@]};do
  fs_array=(${suffix[$rep]})
  for fs in  ${fs_array[@]};do
    java -Xms10G -Xmx20G -jar /mnt/disk4/haitao/software/juicer_tools.2.20.00.jar apa -k NONE -u -r 10000 \
  ${hicdir}/${rep}/${rep}.allValidPairs.hic  ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe \
  ${loopdir}/${rep}_preprocessing_any/perchrAPA/${rep}_${fs} 

    java -Xms10G -Xmx20G -jar /mnt/disk4/haitao/software/juicer_tools.2.20.00.jar apa -k NONE  -u -r 10000 \
  ${hicdir}/HMEC/HMEC.allValidPairs.hic  ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe \
  ${loopdir}/${rep}_preprocessing_any/perchrAPA/HMEC_${fs}
  done;
done
  
  
#loop的peak文件和分开的bedpe文件，即分开的unique，common文件有什么区别吗，可以合并在一起吗后者
#前后的hic和loop文件应该是对应的
#输出文件，${rep}_${fs}，前者是hic，后者是loop文件
#上面代码中原始是无-u的，补充了-u之后call出来的就是perchr，否则是全局

#下面是未交叉脚本
#for rep in ${myrep[@]};do
#  fs_array=(${suffix[$rep]})
#  for fs in  ${fs_array[@]};do
#    java -Xms5G -Xmx10G -jar /mnt/disk4/haitao/software/juicer_tools.2.20.00.jar apa -k NONE -r 10000 \
#  ${hicdir}/${rep}/${rep}.allValidPairs.hic  ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe \
#  ${loopdir}/${rep}_preprocessing_any/${fs}_new 
#  done;
#  for fs in ${HMECsuffix[@]};do
#    java -Xms5G -Xmx10G -jar /mnt/disk4/haitao/software/juicer_tools.2.20.00.jar apa -k NONE  -r 10000 \
#  ${hicdir}/HMEC/HMEC.allValidPairs.hic  ${loopdir}/${rep}_preprocessing_any/loops_${fs}_10000.bedpe \
#  ${loopdir}/${rep}_preprocessing_any/${fs}_new
#  done;
#done