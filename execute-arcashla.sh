#!/bin/bash

# Only specify the R1 fastq files to this script
# Run this script from the ../arcasHLA directory

myfiles=$@ # only specify the R1 fastq files here

echo $myfiles

for myfile in $myfiles; do
  myprefix=`basename $myfile _L001_R1_001.fastq.gz`
  echo $myprefix
  ./arcasHLA genotype /data/volume_2/HLA/${myprefix}_L001_R1_001.fastq.gz /data/volume_2/HLA/${myprefix}_L001_R2_001.fastq.gz -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o /data/volume_2/HLA/output -t 1 -v 
done

echo "FINISHED"
