#merge bam»¹ÊÇbedtools intersect£¬call peak
mamba activate macs2
macs2 callpeak -t SRR13755433.bam SRR13755434.bam -c SRR15567107.bam SRR15567108.bam -f BAM -g hs -n HCC70 -B --outdir HCC70
macs2 callpeak -t SRR13755435.bam SRR13755436.bam -c SRR15567109.bam SRR15567110.bam -f BAM -g hs -n BT549 -B --outdir BT549
macs2 callpeak -t SRR13755425.bam SRR13755426.bam -c SRR15567099.bam SRR15567100.bam -f BAM -g hs -n HMEC -B --outdir HMEC