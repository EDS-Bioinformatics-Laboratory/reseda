#!/bin/bash

celltypes="IGH_HUMAN IGK_HUMAN IGL_HUMAN TRB_HUMAN TRA_HUMAN"

for celltype in $celltypes; do
  python TranslateAndExtractCdr3.py -c $celltype *.assembled.fastq.gz
done

echo "FINISHED"

