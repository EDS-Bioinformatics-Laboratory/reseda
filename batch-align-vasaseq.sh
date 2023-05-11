#!/bin/bash

myfiles=`ls sampled-R2/T*.fastq.gz`
outdir="align-VJC-sampled-R2"
ref="VJC_human.fasta"

mkdir $outdir

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

