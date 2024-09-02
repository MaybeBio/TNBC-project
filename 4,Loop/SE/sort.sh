#!/bin/bash

for file in *.bed;do
  mkdir -p sortbed
  base=$(basename ${file})
  prefix=${base%.*}
  sort -k1,1V -k2,2n -k3,3n ${file} > ${prefix}_sort.bed
  mv ${prefix}_sort.bed sortbed
done