#!/bin/bash

#��ʱ����HMEC�Լ�TNBC�и��Ե�SE��Ϊ�ȶԣ�����Ч��
loopdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result
allEdir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/SE/enhancer/sortbed
Efile=(SE TE allE)
#myrep=(TNBC HMEC) #HMEC BT549 HCC70 MB231 TNBC
declare -A suffix=(
    ["HMEC"]="HMEC_unique HMEC_common HMEC_all" 
    ["TNBC"]="TNBC_unique TNBC_common TNBC_all"
)


#for TNBC,���Ƿ��ļ���ִ�в���
myrep=(TNBC HMEC) #BT549 HCC70 MB231 TNBCͬ��
mkdir -p ${loopdir}/TNBC_preprocessing_any/EPC
for rep in ${myrep[@]};do
  for E in ${Efile[@]};do
    fs_array=(${suffix[$rep]})
    for fs in  ${fs_array[@]};do
      python /mnt/disk4/haitao/bysj_seu/geo_data/hic/script7/mustache_result/EPC/find_functional_loops.py ${allEdir}/${rep}_${E}_sort.bed  ${loopdir}/TNBC_preprocessing_any/loops_${fs}_10000_summit.bed  ${loopdir}/TNBC_preprocessing_any/EPC/${rep}_${E}PC_${fs}.txt
    done;
  done;
done

#������Ը��������suffix��׺��HMEC����������3��TNBC�����ã�Ȼ��rep�Ļ������ø��Եıȶ��У�Ȼ��ע��loop�ļ��Լ�E�ļ���ַ�����ã���������������Է��ڶ�Ӧ���ļ��������EPC