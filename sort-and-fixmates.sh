#!/bin/bash

myfiles=`ls align-VJC-sampled-R2/T*-mapped.fastq.gz`

# Sort fastq files
for f in $myfiles; do
  mydir=`dirname $f`
  prefix=`basename $f .fastq.gz`
  gunzip -c $f |python SortFastq.py |gzip > ${mydir}/${prefix}-sorted.fastq.gz
done

# Fix mates fastq files
myR1=`ls align-VJC-sampled-R2/*_R1_001-mapped-sorted.fastq.gz`
for f in $myR1; do
  mydir=`dirname $f`
  prefix=`basename $f _R1_001-mapped-sorted.fastq.gz`
  python FixMatesFastq.py ${f} ${mydir}/${prefix}_R2_001-mapped-sorted.fastq.gz
done

echo "FINISED"

