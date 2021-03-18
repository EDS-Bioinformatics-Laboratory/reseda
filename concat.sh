ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run17/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20170821-RUN17-datasheet-UMI-DOMINO-new.json -n Domino -c IGH_HUMAN -pre run17-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run17/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20170821-RUN17-datasheet-UMI-DOMINO-new.json -n Domino -c IGH_HUMAN -pre run17-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run17/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20170821-RUN17-datasheet-UMI-DOMINO-new.json -n Domino -c IGH_HUMAN -pre run17-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run19/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20171005_RUN019_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run19-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run19/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20171005_RUN019_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run19-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run19/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20171005_RUN019_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run19-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run20/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20171127_RUN20-Datasheet-UMI-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run20-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run20/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20171127_RUN20-Datasheet-UMI-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run20-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run20/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20171127_RUN20-Datasheet-UMI-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run20-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run21/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20171218_RUN021_Datasheet-umi-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run21-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run21/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20171218_RUN021_Datasheet-umi-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run21-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run21/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20171218_RUN021_Datasheet-umi-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run21-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run22/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180119_RUN22_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run22-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run22/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180119_RUN22_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run22-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run22/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180119_RUN22_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run22-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run23/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180302_RUN23_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run23-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run23/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180302_RUN23_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run23-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run23/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180302_RUN23_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run23-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run26/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180409_RUN26_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run26-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run26/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180409_RUN26_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run26-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run26/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180409_RUN26_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run26-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run28/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180514_RUN28_Datasheet-UMI-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run28-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run28/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180514_RUN28_Datasheet-UMI-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run28-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run28/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180514_RUN28_Datasheet-UMI-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run28-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run29/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180605_RUN29_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run29-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run29/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180605_RUN29_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run29-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run29/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180605_RUN29_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run29-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run30/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180723_RUN30_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run30-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run30/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180723_RUN30_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run30-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run30/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180723_RUN30_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run30-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run31/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180921_RUN31_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run31-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run31/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180921_RUN31_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run31-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run31/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20180921_RUN31_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run31-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run32/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20181125_RUN32_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run32-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run32/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20181125_RUN32_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run32-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run32/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20181125_RUN32_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run32-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run33/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20181214_RUN33_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run33-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run33/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20181214_RUN33_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run33-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run33/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20181214_RUN33_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run33-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run34/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20190214_RUN34_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run34-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run34/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20190214_RUN34_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run34-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run34/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20190214_RUN34_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run34-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run35/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20190609_RUN35_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run35-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run35/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20190609_RUN35_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run35-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run35/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20190609_RUN35_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run35-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run36/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20191003_RUN36_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run36-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run36/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20191003_RUN36_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run36-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run36/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20191003_RUN36_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run36-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run38/final/*IGH_HUMAN*.rr.clones_subs.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20200529_RUN38_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run38-cdr3-clones- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run38/final/*IGH_HUMAN*-all_info.csv.rr.csv > SAMPLES ;  wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20200529_RUN38_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run38-assign-info- $myfiles ; wait
rm $myfiles
#
ls /mnt/immunogenomics/RUNS/project-20210304-DOMINO/processing/20210304-reseda/run38/final/*IGH_HUMAN-clones-mut-sites-reassigned.csv > SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20200529_RUN38_Datasheet-DOMINO-new.json -n DOMINO -c IGH_HUMAN -pre run38-vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
echo 'EVERYTHING DONE'
