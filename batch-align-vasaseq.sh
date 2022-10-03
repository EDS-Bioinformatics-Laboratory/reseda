#!/bin/bash

myfiles=`ls split/*.fastq.gz`
outdir="align-IGLV"
ref="IGLV_human.fasta"

for myfile in $myfiles; do
  prefix=`basename $myfile .fastq.gz`
  echo "bwa mem ${ref} $myfile > ${outdir}/${prefix}.sam"
  bwa mem ${ref} $myfile > ${outdir}/${prefix}.sam
  wait
  samtools view -F 0x04 ${outdir}/${prefix}.sam > align/${prefix}-mapped.sam
  wait
  rm -f ${outdir}/${prefix}.sam
done

echo "FINISHED"

