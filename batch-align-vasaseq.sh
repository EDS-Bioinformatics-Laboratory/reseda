#!/bin/bash

myfiles=`ls split/T*.fastq.gz`

for myfile in $myfiles; do
  prefix=`basename $myfile .fastq.gz`
  echo "bwa mem IGHV_human.fasta $myfile > align/${prefix}.sam"
  bwa mem IGHV_human.fasta $myfile > align/${prefix}.sam
  wait
  samtools view -F 0x04 align/${prefix}.sam > align/${prefix}-mapped.sam
  wait
  rm -f align/${prefix}.sam
done

echo "FINISHED"

