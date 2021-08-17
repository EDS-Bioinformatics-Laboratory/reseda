ls /mnt/immunogenomics/RUNS/project-20210816-AB-ADA-TCRb/processing/20210816-reseda-run28/final/*TRB_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles1=`cat LOCAL_SAMPLES`
ls /mnt/immunogenomics/RUNS/project-20210816-AB-ADA-TCRb/processing/20210816-reseda-run24/final/*TRB_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles2=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20210816-DataSheet-PROJECT-AB-ADA-TCRb.json -n AB-ADA -c TRB_HUMAN -pre cdr3-clones- $myfiles1 $myfiles2 ; wait
rm $myfiles1 $myfiles2


#
ls /mnt/immunogenomics/RUNS/project-20210816-AB-ADA-TCRb/processing/20210816-reseda-run28/final/*TRB_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles1=`cat LOCAL_SAMPLES`
ls /mnt/immunogenomics/RUNS/project-20210816-AB-ADA-TCRb/processing/20210816-reseda-run24/final/*TRB_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles2=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20210816-DataSheet-PROJECT-AB-ADA-TCRb.json -n AB-ADA -c TRB_HUMAN -pre assign-info- $myfiles1 $myfiles2 ; wait
rm $myfiles1 $myfiles2
#


ls /mnt/immunogenomics/RUNS/project-20210816-AB-ADA-TCRb/processing/20210816-reseda-run28/final/*TRB_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles1=`cat LOCAL_SAMPLES`
ls /mnt/immunogenomics/RUNS/project-20210816-AB-ADA-TCRb/processing/20210816-reseda-run24/final/*TRB_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles2=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20210816-DataSheet-PROJECT-AB-ADA-TCRb.json -n AB-ADA -c TRB_HUMAN -pre vjcdr3-clones-mut- $myfiles1 $myfiles2 ; wait
rm $myfiles1 $myfiles2
#
echo 'FINISHED'
