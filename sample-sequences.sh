#!/bin/bash

barcodes=$@
ssize="500 1000 5000 10000 50000 100000 500000 1000000"

# -s seed
for b in $barcodes; do
  for s in $ssize; do
    echo "Barcode ${b} Sample ${s}"

    seqtk sample -s${s} split-R2/${b}_S1_L001_R2_001.fastq.gz ${s} > ${b}-${s}_S1_L001_R2_001.fastq
    wait

    gzip ${b}-${s}_S1_L001_R2_001.fastq
    wait
  done
done

echo "FINISHED"

