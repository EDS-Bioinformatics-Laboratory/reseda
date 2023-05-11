#!/bin/bash

ssize="500 1000 5000 10000 50000 100000 500000 1000000"

# -s seed
for s in $ssize; do
  echo "Sample ${s}"

  seqtk sample -s${s} split-R2/TACTGGTA_S1_L001_R2_001.fastq.gz ${s} > TACTGGTA-${s}_S1_L001_R2_001.fastq
  wait

  gzip TACTGGTA-${s}_S1_L001_R2_001.fastq
  wait
done

echo "FINISHED"

