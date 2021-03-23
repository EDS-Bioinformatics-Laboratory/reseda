#!/bin/bash

myfiles=`ls SAMPLES-run*`
ref=reference.fasta

for myfile in $myfiles; do
    # Retrieve the fastq files
    cp $myfile SAMPLES
    ./copy-from-beehub.sh
    wait

    # Perform pairwise assembly
    r1_samples=`grep R1_001 LOCAL_SAMPLES`
    ./batch-pear.sh ${r1_samples}
    wait

    # Perform the alignment
    samples=`ls *.assembled.fastq.gz`
    ./batch-align.sh ${ref} ${samples}
    wait

    # Transfer to the researchdrive
    myurl=https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/PROJECTS/20210218-dornatien-SHM-mouse/processing/20210301-alignment/Results/
    ./copy-to-webdav.sh ${myurl}/Alignments/ *.sam
    wait
    ./copy-to-webdav.sh ${myurl}/SNPs/ *.snp.csv
    wait

    # Remove fastq and alignment files before the next batch
    rm -f *L001*
    wait
done

echo "FINISHED"

