#!/bin/bash

mydir="align-IGLV"
myref="IGLV_human.fasta"
myfiles=`ls ${mydir}/*.sam`

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

