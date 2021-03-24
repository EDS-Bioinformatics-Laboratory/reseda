#
ls /mnt/immunogenomics/RUNS/run35-20190609-miseq/no-umi/final/*IGH_MOUSE*-all_info.csv.rr.all_info.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20190609_RUN35_Datasheet-VDJmouse-new.json -n VDJmouse -c IGH_MOUSE -pre all-info- $myfiles ; wait
rm $myfiles
echo 'FINISHED'
