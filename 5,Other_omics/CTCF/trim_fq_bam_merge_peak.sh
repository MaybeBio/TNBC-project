#!/bin/bash
mamba activate qc
fqdir=/mnt/disk4/haitao/bysj_seu/geo_data/other_omics/CTCF/rawfq
fqs=(SRR13755426 SRR13755436 SRR15567100 SRR15567110 SRR13755425 SRR13755434 SRR15567099 SRR15567108 SRR13755433 SRR13755435 SRR15567107 SRR15567109) #12

#fq ����,trim1
for fq in ${fqs[@]};do
  trim_galore --fastqc_args "-o multiqc1 -t 15 " -j 15 ${fqdir}/${fq}.fastq
  mv ${fq}_trimmed.fq  ${fq}.fastq
done
multiqc ./multiqc1

#trim2,���ݾ���ǰ���fqӦ����3'�ζ�2p�������һ�����Ҫ����ǰ�󱨸�������
for fq in ${fqs[@]};do
  trim_galore --fastqc_args "-o multiqc2 -t 15" --three_prime_clip_R1 2 -j 15  ${fq}.fastq 
  mv ${fq}_trimmed.fq  ${fq}.fq
done
multiqc ./multiqc2

#ʣ�µ�${fq}_trimmed.fq�ļ��������յ�call bam,bam�����������õ���${fq}.bam�ļ�
for fq in ${fqs[@]};do
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
 
#merge bam����bedtools intersect��call peak
mamba activate macs2
macs2 callpeak -t SRR13755433.bam SRR13755434.bam -c SRR15567107.bam SRR15567108.bam -f BAM -g hs -n HCC70 -B --outdir HCC70
macs2 callpeak -t SRR13755435.bam SRR13755436.bam -c SRR15567109.bam SRR15567110.bam -f BAM -g hs -n BT549 -B --outdir BT549
macs2 callpeak -t SRR13755425.bam SRR13755426.bam -c SRR15567099.bam SRR15567100.bam -f BAM -g hs -n HMEC -B --outdir HMEC