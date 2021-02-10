#!/bin/bash

uploadurl="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/runXXX-20210210-briney/data/raw/"

urls=`cat _LINKS.md`

for url in $urls; do
    myfile=`basename $url`
    echo $url
    echo $myfile

    wget $url
    wait

    ./copy-to-webdav.sh $uploadurl $myfile
    wait

    rm $myfile
    wait
done

echo "FINISHED"
