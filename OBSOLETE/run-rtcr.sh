#!/bin/bash

starttime=`date +%s`

samples=$@  # fastq files

thisdir=`pwd`
for s in $samples; do
    prefix=`basename $s .fastq.gz`
    mkdir $prefix
    cd $prefix
    rtcr run -t 2 -i ../$s
    cd $thisdir
    wait
done

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"
