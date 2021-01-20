#!/bin/bash

myfiles=`cat SAMPLES-all-contig-annotations`

for myfile in $myfiles; do
    echo $myfile > SAMPLES
    ./copy-from-beehub.sh
    wait

    samplename=`echo $myfile | perl -ne '@c=split(/\//); print "$c[-3]\n"'`
    mv all_contig_annotations.csv $samplename-all_contig_annotations.csv
    echo $samplename-all_contig_annotations.csv
done

echo "FINISHED"

