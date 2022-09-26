#!/bin/bash

starttime=`date +%s`

samples=$@  # get all arguments

# TTTGTGTC_S1_L001_R1_001-mapped.fastq.gz

for s in $samples; do
    mydir=`dirname $s`
    prefix=`basename $s _L001_R1_001-mapped.fastq.gz`
    ./pear-0.9.6-bin-64 -f ${mydir}/${prefix}_L001_R1_001-mapped.fastq.gz -r ${mydir}/${prefix}_L001_R2_001-mapped.fastq.gz -o ${prefix} > ${prefix}-pear.log 2> ${prefix}-pear.err
    gzip -f ${prefix}.assembled.fastq
done

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"
