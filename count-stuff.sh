#!/bin/bash

barcodes="AAGCGAGT CACCATGT CAGCTTTG CCGTCTAA CTTACGAC GAAACTGC GGGGATTT TACCAGAC TACTGGTA TCCGAACA TCCTAGCT"

for b in $barcodes; do
  wc -l align-VJC-sampled-R2/${b}-*-mapped.sam > counts-sampled/${b}-count-alignments.txt
  grep "Assembled reads \." assembly-sampled/${b}-*-pear.log > counts-sampled/${b}-count-assembled.txt
  grep '^4' cdr3-sampled/${b}-*-report.txt > counts-sampled/${b}-count-cdr3.txt
done

echo "FINISHED"

