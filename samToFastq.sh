#!/bin/bash

mydir="align-VJC-sampled-R2"
myref="VJC_human.fasta"
myfiles=`ls ${mydir}/*-mapped.sam`

for sam in ${myfiles}; do
  prefix=`basename ${sam} .sam`
  echo ${prefix}

  samtools view -bh -T ${myref} -o ${mydir}/${prefix}.bam ${sam}
  wait

  samtools fastq ${mydir}/${prefix}.bam > ${mydir}/${prefix}.fastq
  wait

  gzip ${mydir}/${prefix}.fastq
  wait

  rm -f ${mydir}/${prefix}.bam
done

echo "FINISHED"

