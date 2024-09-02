#!/bin/bash
fqdir=/mnt/disk4/haitao/bysj_seu/geo_data/other_omics/CTCF/rawfq
fqs=(SRR13755426 SRR13755436 SRR15567100 SRR15567110)

#下面是第一次trim
#for fq in ${fqs[@]};do
#  trim_galore --fastqc_args "-o multiqc -t 10 " -j 15 ${fq}.fastq
#done


for fq in ${fqdir}/${fqs[@]};do
  trim_galore --fastqc_args "-o multiqc2 -t 15" --three_prime_clip_R1 2 -j 15  ${fq}_trimmed.fq 
done
multiqc ./multiqc2