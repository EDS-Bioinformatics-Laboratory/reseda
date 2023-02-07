#!/bin/bash

mydir=assembly-VJC-R2

python TranslateAndExtractCdr3.py -c IGH_HUMAN ${mydir}/*.fastq.gz
wait

python TranslateAndExtractCdr3.py -c IGK_HUMAN ${mydir}/*.fastq.gz
wait

python TranslateAndExtractCdr3.py -c IGL_HUMAN ${mydir}/*.fastq.gz
wait

python TranslateAndExtractCdr3.py -c TRB_HUMAN ${mydir}/*.fastq.gz
wait

python TranslateAndExtractCdr3.py -c TRA_HUMAN ${mydir}/*.fastq.gz
wait

echo "FINISHED"

