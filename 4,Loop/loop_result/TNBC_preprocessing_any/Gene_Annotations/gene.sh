#!/bin/bash

for file in *.bed;do
  base=$(basename ${file})
  prefix=${base%.*}
  cut -f 1 ${file} > ${prefix}_GSEA.txt
done