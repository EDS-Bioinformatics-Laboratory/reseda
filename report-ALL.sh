#!/bin/bash

run="run10-20161031-miseq"
resultsdir="results-tbcell-nov2016"

files="/mnt/immunogenomics/RUNS/${run}/${resultsdir}/reports/*.log"
grep '^Assembled reads \.' $files > report-PEAR.txt
cat report-PEAR.txt

files="/mnt/immunogenomics/RUNS/${run}/${resultsdir}/reports/*-midcount.txt"
more $files > report-MIDs.txt
cat report-MIDs.txt

files="/mnt/immunogenomics/RUNS/${run}/${resultsdir}/reports/*HUMAN-report.txt /mnt/immunogenomics/RUNS/${run}/${resultsdir}/reports/*MOUSE-report.txt"
grep '^4' $files > report-CDR3.txt
cat report-CDR3.txt

files="/mnt/immunogenomics/RUNS/${run}/${resultsdir}/reports/*-productive.txt"
grep 'Unique reads with V and J in all_info' $files > report-PRODUCTIVE.txt
cat report-PRODUCTIVE.txt

python report-after-v-reassignment.py /mnt/immunogenomics/RUNS/${run}/${resultsdir}/final/correct-mid/*.rr.clones_subs.csv
cat report-AFTER-V-REASSIGNMENT.txt

python report-combine-all.py
wait

echo "FINISHED"
