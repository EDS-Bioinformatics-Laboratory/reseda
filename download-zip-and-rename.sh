#!/bin/bash

urls=`cat SAMPLES-cellranger-output`

for url in $urls; do
    echo $url > SAMPLES

    sample=`echo $url |perl -ne '@c=split(/\//); print "$c[-3]\n";'`
    echo $sample    

    ./copy-from-beehub.sh
    wait

    newsample="${sample}_S0_L001.assembled.fastq"
    mv all_contig.fastq $newsample
    wait

    gzip $newsample
    wait
done

echo "FINISHED"

