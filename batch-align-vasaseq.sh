#!/bin/bash

myfiles=`ls *.fastq.gz`
outdir="align-ALL-V"
ref="ALL_V_human.fasta"

for myfile in $myfiles; do
  prefix=`basename $myfile .fastq.gz`
  echo "bwa mem ${ref} $myfile > ${outdir}/${prefix}.sam"
  bwa mem ${ref} $myfile > ${outdir}/${prefix}.sam
  wait
  samtools view -F 0x04 ${outdir}/${prefix}.sam > ${outdir}/${prefix}-mapped.sam
  wait
  rm -f ${outdir}/${prefix}.sam
done

echo "FINISHED"

