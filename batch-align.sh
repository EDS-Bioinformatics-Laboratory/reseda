#!/bin/bash

starttime=`date +%s`

ref=$1; shift  # e.g. IGHV-human.fasta
samples=$@  # get rest of the arguments

for s in $samples; do
#    prefix=`basename $s _R1_001.fastq.gz`
#    ./align-pairs.sh $ref ${prefix}_R1_001.fastq.gz ${prefix}_R2_001.fastq.gz
    ./align-sequences.sh $ref $s
done

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"
