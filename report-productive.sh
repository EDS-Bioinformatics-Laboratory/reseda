files="/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/reports/*-productive.txt"

grep 'Unique reads with V and J in all_info' $files > report-PRODUCTIVE.txt
cat report-PRODUCTIVE.txt
