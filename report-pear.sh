files="/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/reports/*.log"

grep '^Assembled reads \.' $files > report-PEAR.txt
cat report-PEAR.txt
