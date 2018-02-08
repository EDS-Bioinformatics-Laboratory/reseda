#!/bin/bash

# Example URL: myurl="https://beehub.nl/amc-immunogenomics/RUNS/run246/data/"
# Example URL: myurl="https://beehub.nl/amc-immunogenomics/RUNS/run246/results-tbcell/"
# Example URL: myurl="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/run246/data/"
myurl=$1; shift

myfiles="*.rr.* *mutations*"
filelist=`echo $myfiles | perl -ne "@c=split(/\s/); print join(',', @c);"`

starttime=`date +%s`

echo "Calculating checksums..."
shasum $myfiles > CHECKSUM.$starttime
wait

echo "Transferring files..."
curl -T "{$filelist,CHECKSUM.$starttime}" --netrc $myurl
wait

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`

echo "Finished in $difftime seconds"

# Transfer files FROM beehub (-C - is for resuming files)
# curl -C - -L -u bschaik -O https://beehub.nl/home/bschaik/immunology/RUNS/run234/data/blah.sff
