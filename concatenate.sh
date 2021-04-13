ls /mnt/immunogenomics/RUNS/run40-20201116-miseq/20210407-reseda-race-lightchains/final/*IGK_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20201116-RUN40-DataSheet-UMI-Constant-RACE-IGK-new.json -n RACE -c IGK_HUMAN -pre cdr3-clones-IGK- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/run40-20201116-miseq/20210407-reseda-race-lightchains/final/*IGK_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20201116-RUN40-DataSheet-UMI-Constant-RACE-IGK-new.json -n RACE -c IGK_HUMAN -pre assign-info-IGK- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/run40-20201116-miseq/20210407-reseda-race-lightchains/final/*IGK_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20201116-RUN40-DataSheet-UMI-Constant-RACE-IGK-new.json -n RACE -c IGK_HUMAN -pre vjcdr3-clones-mut-IGK- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/run40-20201116-miseq/20210407-reseda-race-lightchains/final/*IGL_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20201116-RUN40-DataSheet-UMI-Constant-RACE-IGL-new.json -n RACE -c IGL_HUMAN -pre cdr3-clones-IGL- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/run40-20201116-miseq/20210407-reseda-race-lightchains/final/*IGL_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20201116-RUN40-DataSheet-UMI-Constant-RACE-IGL-new.json -n RACE -c IGL_HUMAN -pre assign-info-IGL- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/run40-20201116-miseq/20210407-reseda-race-lightchains/final/*IGL_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20201116-RUN40-DataSheet-UMI-Constant-RACE-IGL-new.json -n RACE -c IGL_HUMAN -pre vjcdr3-clones-mut-IGL- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
