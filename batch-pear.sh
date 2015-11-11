#!/bin/bash

starttime=`date +%s`

samples=$@  # get all arguments

for s in $samples; do
    prefix=`basename $s _R1_001.fastq.gz`
    ./pear-0.9.6-bin-64 -f ${prefix}_R1_001.fastq.gz -r ${prefix}_R2_001.fastq.gz -o ${prefix} > ${prefix}-pear.log 2> ${prefix}-pear.err
done

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"
