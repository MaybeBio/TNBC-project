#!/bin/bash

diffTADfolder=/mnt/disk4/haitao/bysj_seu/geo_data/hic/script4/diffTAD
myres=(40000 20000 10000)
myrep=(BT549 HCC70 MB231)
for res in ${myres[@]};do
  for rep in ${myrep[@]};do 
     wc -l ${res}_${rep}_diff/*_tad.bed >> diffTAD.txt
     wc -l ${res}_${rep}_diff/test_to_control_* >> diffTAD.txt
     echo -e 'separate' >> diffTAD.txt
     cut -d ',' -f 1 ${res}_${rep}_diff/test_to_control_separate.csv | sort | uniq -c | wc -l >> diffTAD.txt
     echo -e "${res}_${rep}\n" >> diffTAD.txt
  done;
done
  