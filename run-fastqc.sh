#!/bin/bash

starttime=`date +%s`

# Input
files=$@

# Outdir
outdir="fastqc-reports"

# Run FastQC
./FastQC/fastqc --extract $files

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH FASTQC IN $difftime seconds"
