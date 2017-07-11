#!/bin/bash

cdr3_files=`ls final/correct-mid/*IGH*.rr.all_info.csv`

for cdr3_file in ${cdr3_files}; do
    SAMPLEMID=`basename ${cdr3_file} |perl -ne '@c=split(/-IGH_/);print "$c[0]";'`
    echo ${SAMPLEMID}
    v_file=`ls ${SAMPLEMID}-IGHV*.sam.mut.txt`
    j_file=`ls ${SAMPLEMID}-IGHJ*.sam.mut.txt`

    Rscript --vanilla MutationAnalysisVJ.R indir="." outdir="final/correct-mid" V.file="${v_file}" J.file="${j_file}" CDR3.file="${cdr3_file}"
done
