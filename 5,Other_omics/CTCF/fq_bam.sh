#!/bin/bash
fqfile=(SRR13755426 SRR13755436 SRR15567100 SRR15567110)
for fq in ${fqfile[@]};do
  bowtie2 -X 100 --very-sensitive -p 15 -x /mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19 -U ${fq}.fq | \
  samtools view -bS -@ 15 | \
  samtools sort -@ 15  -O BAM -o ${fq}_pre.bam
  java -Xmx4g -jar ~/software/picard.jar  MarkDuplicates \
    I=${fq}_pre.bam \
    O=${fq}.bam \
    REMOVE_DUPLICATES=true \
    M=${fq}_marked_dup_metrics.txt
  samtools index -@ 15 ${fq}.bam 
done
