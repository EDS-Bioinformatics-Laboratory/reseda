#!/bin/bash

# Example URL: myurl="https://beehub.nl/amc-immunogenomics/RUNS/run246/data/"
# Example URL: myurl="https://beehub.nl/amc-immunogenomics/RUNS/run246/results-tbcell/"
myurl="https://beehub.nl/amc-immunogenomics/RUNS/runNN-2016MMDD-miseq/results-tbcell/reports/"

#myfiles="runmetrics* *.sff"
myfiles="*-pear.log *-pear.err"
filelist=`echo $myfiles | perl -ne "@c=split(/\s/); print join(',', @c);"`

starttime=`date +%s`

echo "Calculating checksums... (you need to type the password for Beehub later, so please wait for approx 30 seconds) ..."
shasum $myfiles > CHECKSUM.$starttime
wait

echo "Transferring files..."
curl -T "{$filelist,CHECKSUM.$starttime}" -u bschaik $myurl
wait

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`

echo "Finished in $difftime seconds"

# Transfer files FROM beehub (-C - is for resuming files)
# curl -C - -L -u bschaik -O https://beehub.nl/home/bschaik/immunology/RUNS/run234/data/blah.sff
