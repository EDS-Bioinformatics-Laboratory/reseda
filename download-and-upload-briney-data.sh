#!/bin/bash

uploadurl="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/runXXX-20210210-briney/data/raw/"
$scriptdir=`pwd`

urls=`cat _LINKS.md`

cd /data2

for url in $urls; do
    myfile=`basename $url`
    echo $url
    echo $myfile

    wget $url
    wait

    tar zxvf $myfile
    wait

    find data -type f -regex '.*.fastq.gz' -exec mv {} . \;
    wait

    $scriptdir/copy-to-webdav.sh $uploadurl *.fastq.gz
    wait

    rm -rf /data2/*
    wait
done

echo "FINISHED"
