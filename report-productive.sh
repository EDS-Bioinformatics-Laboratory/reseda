files="/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/*-productive.txt"

grep 'Productive' $files > report-PRODUCTIVE.txt
cat report-PRODUCTIVE.txt
