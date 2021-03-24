#
ls /mnt/immunogenomics/RUNS/run40-20201116-miseq/reseda20201116/final/*IGH_MOUSE*-all_info.csv.rr.all_info.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20201116-RUN40-DataSheet-VDJmouse-new.json -n VDJmouse -c IGH_MOUSE -pre all-info- $myfiles ; wait
rm $myfiles
echo 'FINISHED'
