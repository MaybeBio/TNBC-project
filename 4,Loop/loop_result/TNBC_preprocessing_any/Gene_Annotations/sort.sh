#!/bin/bash

for file in *.bed;do
  base=$(basename ${file})
  prefix=${base%.*}
  sort -k2,2nr ${file} >  ${prefix}_sort.bed
done