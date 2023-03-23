#!/bin/bash

# align-VJC-R2/TTTTCGCC_S1_L001_R1_001-mapped-ufis.csv
myfiles=`ls align-VJC-R2/*-mapped-ufis.csv`

for f in ${myfiles}; do
  barcode=`basename ${f} _S1_L001_R1_001-mapped-ufis.csv`
  ufis=`cut -f 2 ${f}|sort|uniq |wc -l`
  echo "${barcode} ${ufis}"
done

