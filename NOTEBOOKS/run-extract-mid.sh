#!/bin/bash

samples=`cat SAMPLES-fastq`

for sample in $samples; do
    echo $sample
    echo $sample > SAMPLES
    ../copy-from-beehub.sh
    wait

    local_sample=`cat LOCAL_SAMPLES`
    python2 ../ExtractFirstNucleotides.py $local_sample > ${local_sample}-mids.csv
    wait

    echo "Removing ${local_sample}"
    rm $local_sample
done

echo "FINISHED"

