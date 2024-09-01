#!/bin/bash

# mamba activate RNA
fq1=$1
# fq2=$2
rep_folder=$2

#如果是双端测序，就选用下列命令
#salmon quant --gcBias -l A -1 ${fq1} -2 ${fq2}  -i /mnt/disk4/haitao/bysj_seu/ref_genome/ucsc_hg19/hg19.salmon_sa_index -g /mnt/disk4/haitao/bysj_seu/ref_genome/ucsc_hg19/hg19.gtf -o ${rep_folder} -p 15

#如果是单端测序，就选用下列命令：
# salmon quant --gcBias -l A -r ${fq1}   -i /mnt/disk4/haitao/bysj_seu/ref_genome/ucsc_hg19/hg19.salmon_sa_index -g /mnt/disk4/haitao/bysj_seu/ref_genome/ucsc_hg19/hg19.gtf -o ${rep_folder} -p 15


#下面是针对ucsc的fa+ensemble的gtf

salmon quant --gcBias -l A -r ${fq1}   -i /mnt/disk4/haitao/bysj_seu/ref_genome/ucsc_hg19/hg19.salmon_sa_index -g /mnt/disk4/haitao/bysj_seu/ref_genome/ucsc_hg19/hg19.gtf -o ${rep_folder} -p 20