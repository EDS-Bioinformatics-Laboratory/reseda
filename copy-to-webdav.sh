#!/bin/bash

# Example URL: myurl="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/run246/Data/NameOfDataset_1/Raw/"

myurl=$1; shift
myfiles=$@

if [ "${myfiles}" = "" ]; then
    echo "Usage: ./copy-to-webdav.sh WEBDAVURL FILE1 [FILE2]"
    exit
fi

filelist=`echo $myfiles | perl -ne "@c=split(/\s/); print join(',', @c);"`

starttime=`date +%s`

echo "Calculating checksums..."
shasum $myfiles > CHECKSUM.$starttime
wait

echo "Transferring files to ${myurl} ..."
curl -T "{$filelist,CHECKSUM.$starttime}" --netrc $myurl
wait

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`

echo "Finished in $difftime seconds"

# Transfer files FROM beehub (-C - is for resuming files)
# curl -C - -L -u bschaik -O https://beehub.nl/home/bschaik/immunology/RUNS/run234/data/blah.sff
