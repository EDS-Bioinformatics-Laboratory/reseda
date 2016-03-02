files="/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/reports/*HUMAN-report.txt"

grep '^4' $files > report-CDR3.txt
cat report-CDR3.txt
