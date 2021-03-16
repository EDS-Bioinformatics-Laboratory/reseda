#!/bin/bash

myfiles=`ls *sheet*.csv`

for myfile in $myfiles; do
    jsonfile=`basename $myfile .csv`.json
    echo $myfile $jsonfile
    python MetaData.py $myfile > $jsonfile
done

echo "FINISHED"

