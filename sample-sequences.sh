#!/bin/bash

ssize="500 1000 5000 10000 50000 100000 500000 1000000"

# -s seed
for s in $ssize; do
  echo "Sample ${s}"

  seqtk sample -s${s} split-R1/AAGCGAGT_S1_L001_R1_001.fastq.gz ${s} > AAGCGAGT-${s}_S1_L001_R1_001.fastq
  seqtk sample -s${s} split-R2/AAGCGAGT_S1_L001_R2_001.fastq.gz ${s} > AAGCGAGT-${s}_S1_L001_R2_001.fastq
  wait

  gzip AAGCGAGT-${s}_S1_L001_R1_001.fastq
  gzip AAGCGAGT-${s}_S1_L001_R2_001.fastq
  wait
done

echo "FINISHED"

