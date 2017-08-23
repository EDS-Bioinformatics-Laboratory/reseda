#!/bin/bash

run="run16-20170808-miseq"
runinfo="20170804_RUN16_Datasheet-new.json"
resultsdir="results-tbcell"

files="/mnt/immunogenomics/RUNS/${run}/${resultsdir}/reports/*.log"
grep '^Assembled reads \.' $files > report-PEAR.txt
echo "Wrote report-PEAR.txt"

files="/mnt/immunogenomics/RUNS/${run}/${resultsdir}/reports/*-midcount.txt"
more $files > report-MIDs.txt
echo "Wrote report-MIDs.txt"

files="/mnt/immunogenomics/RUNS/${run}/${resultsdir}/reports/*HUMAN-report.txt /mnt/immunogenomics/RUNS/${run}/${resultsdir}/reports/*MOUSE-report.txt"
grep '^4' $files > report-CDR3.txt
echo "Wrote report-CDR3.txt"

files="/mnt/immunogenomics/RUNS/${run}/${resultsdir}/reports/*-productive.txt"
grep 'Unique reads with V and J in all_info' $files > report-PRODUCTIVE.txt
echo "Wrote report-PRODUCTIVE.txt"

files="/mnt/immunogenomics/RUNS/${run}/${resultsdir}/raw/*.sam"
python report-alignments.py ${files}
echo "Wrote report-ALIGNED-*"

python2 report-after-v-reassignment.py /mnt/immunogenomics/RUNS/${run}/${resultsdir}/final/correct-mid/*.rr.clones_subs.csv
echo "Wrote report-AFTER-V-REASSIGNMENT.txt"

python2 report-combine-all.py ${runinfo}
wait

echo "FINISHED"
