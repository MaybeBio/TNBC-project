import sys
from os import system
name=sys.argv[1]
system(f"macs2 callpeak -t {name}.bam -f BAMPE -g hs -n {name} --outdir {name} -B --shift 100 --extsize 200 -p 0.001 --seed 2024")
