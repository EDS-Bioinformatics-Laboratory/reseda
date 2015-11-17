files="/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/reports/*.log"

grep '^Assembled reads \.' $files > report-PEAR.txt
cat report-PEAR.txt
