#!/bin/bash

# Preparation:
# mv reference/* .
# mv hlaforest/scripts/* .

# Create a file with all the fastq files:
# ls TESTDATA/* > SAMPLES

# Configure this:
celltype="IGH_HUMAN"
mids="MIDS-miseq.txt"
refs="IGHV_human.fasta IGHJ_human.fasta"
v="IGHV_human"
j="IGHJ_human"
beehub_mount="/mnt/immunogenomics/RUNS/barbera-test"
beehub_web="https://beehub.nl/amc-immunogenomics/RUNS/barbera-test"

# Then run ./execute-all.sh


#######################

starttime=`date +%s`
ip_address=`hostname -I`

function test {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "ERROR with $1" >&2
        cat SAMPLES | mail -s "ERROR ${ip_address} ${celltype}" b.d.vanschaik@amc.uva.nl
        exit
    fi
    return $status
}

samples=`cat SAMPLES`  # get all arguments
r1_samples=`grep R1_001 SAMPLES`

### Analysis on raw fastq files ###

# FastQC
test ./run-fastqc.sh ${samples}
wait


# Pairwise assembly
test ./batch-pear.sh ${r1_samples}
wait

### Continue with assembled fastq files ###

samples=`ls *.assembled.fastq.gz`

# Split on MID
test python fastq-split-on-mid.py ${mids} split ${samples}
wait

### Continue with the assembled, split per mid, fastq files ###

samples=`ls split/*.fastq.gz`

# FastQC report
test ./run-fastqc.sh ${samples}
wait

# Search for primers in the fastq files
test python motif-search-batch.py ${samples}
wait

# Extract the CDR3 sequence
test python translate-and-extract-cdr3.py ${celltype} ${samples}
wait

# Align sequences against IMGT and call SNPs
for ref in $refs; do
    test ./batch-align.sh ${ref} ${samples}
done
wait

### Continue with the aligned sequences ###

bamfiles=`ls *.sam`

# Alignment quality report TO IMPLEMENT

# HLAforest for HLA samples
# for bam in $bamfiles; do
#     ./hla-forest.sh ${ref} ${bam}
#     wait
# done

### Generate reports ###

mkdir final

# For each sample; do
for sample in ${samples}; do
    mydir=`dirname ${sample}`
    prefix=`basename ${sample} .fastq.gz`

    # Combine MID, CDR3, V, J and sequence information
    midFile=`echo ${mydir}/${prefix}|perl -ne 's/(.+)-.+$/$1-report.txt/;print;'`
    cdr3File=${sample}-${celltype}-CDR3.csv
    vFile=${prefix}-${v}-easy-import.txt
    jFile=${prefix}-${j}-easy-import.txt
    seqFile=${sample}-${celltype}.csv
    outFile="final/${prefix}-${celltype}-all_info.csv"
    cloneFile="final/${prefix}-${celltype}-clones.csv"
    totalFile="final/${prefix}-${celltype}-productive.txt"
    test python combine-immuno-data.py ${midFile} ${cdr3File} ${vFile} ${jFile} ${seqFile} ${outFile} ${cloneFile} ${totalFile}
    wait
done

# Make output directories
mkdir ${beehub_mount}/results-tbcell
mkdir ${beehub_mount}/results-tbcell/raw
mkdir ${beehub_mount}/results-tbcell/reports
mkdir ${beehub_mount}/results-tbcell/final

# Transfer data to Beehub
test ./copy-to-beehub-reports.sh ${beehub_web}/results-tbcell/reports/
test ./copy-to-beehub-raw.sh ${beehub_web}/results-tbcell/raw/
cd split
test ./copy-to-beehub-reports.sh ${beehub_web}/results-tbcell/reports/
test ./copy-to-beehub-raw.sh ${beehub_web}/results-tbcell/raw/
cd ../final
test ./copy-to-beehub-reports.sh ${beehub_web}/results-tbcell/reports/
test ./copy-to-beehub-final.sh ${beehub_web}/results-tbcell/final/
cd ..

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"

# Send mail to Barbera to tell that the analysis is finished
cat SAMPLES | mail -s "FINISHED ${ip_address} ${celltype} - ${difftime} seconds" b.d.vanschaik@amc.uva.nl
