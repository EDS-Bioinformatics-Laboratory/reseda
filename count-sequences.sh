#!/bin/bash

myfiles=`ls split-R2/*.fastq.gz`

for x in $myfiles; do
  barcode=`basename $x _S1_L001_R2_001.fastq.gz`
  count=`gunzip -c $x|grep '@'|wc -l`
  echo "$barcode $count"
done

