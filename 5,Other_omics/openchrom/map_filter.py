import sys  
from os import system  

fq1,fq2=sys.argv[1],sys.argv[2]
name=sys.argv[3]

def run_map_filter_callpeak(fq1,fq2,name):
    system(f"bowtie2 -X 100 --very-sensitive -p 10 -1 {fq1} -2 {fq2} -x /mnt/disk4/haitao/bysj_seu/ref_genome/new_hg19/hg19 | \
    samtools view -bS -@ 10 | \
    samtools sort -@ 10  -O BAM -o {name}_pre.bam ")
    system(f"java -Xmx4g -jar ~/software/picard.jar  MarkDuplicates \
    I={name}_pre.bam \
    O={name}.bam \
    REMOVE_DUPLICATES=true \
    M={name}_marked_dup_metrics.txt ")
    system(f"samtools index -@ 10 {name}.bam ")

run_map_filter_callpeak(fq1,fq2,name)


  
    
