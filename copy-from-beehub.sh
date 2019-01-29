#!/bin/bash

basemnt="\/mnt\/immunogenomics\/"
baseurl="https:\/\/researchdrive.surfsara.nl\/remote.php\/webdav\/amc-immunogenomics\/"

starttime=`date +%s`

samples=`cat SAMPLES`  # get all samples: /mnt/immunogenomics/RUNS/runXXX/data/*.fastq.gz

echo "Transferring files..."
rm -f LOCAL_SAMPLES
for sample in $samples; do
    sampleurl=`echo ${sample} | perl -ne "s/$basemnt/$baseurl/;print;"`
    echo ${sampleurl}
    localsample=`basename ${sample}`
    curl -C - -L --netrc -O $sampleurl
    echo ${localsample} >> LOCAL_SAMPLES
    wait
done

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`

echo "Finished in $difftime seconds"

# Transfer files FROM beehub (-C - is for resuming files)
# curl -C - -L -u bschaik -O https://beehub.nl/home/bschaik/immunology/RUNS/run234/data/blah.sff
