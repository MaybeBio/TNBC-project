#!/bin/bash

#mamba activate hicexplorer

myrep=(HMEC BT549 HCC70 MB231)
mcooldir=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script5/merge_fq_now #mcool file
myres=(40000 20000)       #��GENOVA�е�PESCAn��res�����󣬽���ʹ�õ�res��20kb-40kb����,PE-SCAn 20kb-40kb
#�����������������϶�����HMECΪ�� 
bed1=$1 #����Ƚϵĵ�һ�������bed�ļ�
bed2=$2 #����Ƚϵĵڶ��������bed�ļ�
range=$3

for res in ${myres[@]};do
  for rep in ${myrep[@]};do
  
  base1=$(basename ${bed1})
  base2=$(basename ${bed2})
  prefix_bed1=${base1%.*}
  prefix_bed2=${base2%.*}
  
    hicAggregateContacts  -m ${mcooldir}/${rep}.mcool::/resolutions/${res} --outFileName ${res}_${rep}_${prefix_bed1}_${prefix_bed2}_${range}_none.pdf --BED ${bed1} --BED2 ${bed2}  --operationType mean --mode intra-chr --largeRegionsOperation center --range ${range} 
    hicAggregateContacts  -m ${mcooldir}/${rep}.mcool::/resolutions/${res} --outFileName ${res}_${rep}_${prefix_bed1}_${prefix_bed2}_${range}.pdf --BED ${bed1} --BED2 ${bed2} --transform obs/exp --operationType mean --mode intra-chr --largeRegionsOperation center  --range ${range} 
    
    done;
done
    
#--considerStrandDirection ���������CTCF��ʱ����ʹ��
#--outFilePrefixMatrix  --outFileContactPairs --outFileObsExp --diagnosticHeatmapFile �������Щѡ���ȫ����ʱ��ʹ�ã�����������������
#��֮ʹ������bed�ļ��Ŀ��Ժ�GENOVA�е�CSCAn���߻�����Ȼ��ʹ��һ��bed�ļ��Ŀ��Ժ�GENOVA�е�PESCAn���߻���
#--range ��ʱѡ��200kb-20mb