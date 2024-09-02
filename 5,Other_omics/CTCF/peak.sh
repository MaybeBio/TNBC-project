#!/bin/bash

macs2 callpeak -t SRR13755436.bam -c SRR15567110.bam -f BAM -g hs -n BT549_rep2 -B  

macs2 callpeak -t SRR13755426.bam -c SRR15567100.bam -f BAM -g hs -n HMEC_rep2 -B

#貌似macs2 callpeak可以接受多个输入，不知道能不能直接这里就按照rep1+rep2来操作