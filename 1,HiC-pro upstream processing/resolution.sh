#!/bin/bash

for i in *.allValidPairs
do
prefix=$(basename $i .allValidPairs)

echo "$prefix"
bash  /mnt/disk4/haitao/software/juicer-1.6/misc/calculate_map_resolution.sh $i ${prefix}.txt

done