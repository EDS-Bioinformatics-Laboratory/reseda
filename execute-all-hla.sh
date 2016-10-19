#!/bin/bash

# Preparation:
# mv reference/* .
# mv hlaforest/scripts/* .

# Create a file with all the fastq files:
# ls TESTDATA/* > SAMPLES

# Get settings from commandline arguments
run=$1
mids=$2
organism=$3
cell=$4
celltype=$5

# Configure this:
# cell="HLA"
# organism="human"
# celltype="${cell}_HUMAN"
# run="run09-20160919-miseq"

# Reference sequences
#mids="MIDS-miseq.txt"
#refs="${cell}V_${organism}.fasta ${cell}J_${organism}.fasta"
refs="hla_nuc_nospace.fasta"
v="${cell}V_${organism}"
j="${cell}J_${organism}"

# Mount the Beehub webdav server and configure the location
beehub_mount="/mnt/immunogenomics/RUNS/${run}"
beehub_web="https://beehub.nl/amc-immunogenomics/RUNS/${run}"

# Then run ./execute-all.sh


#######################

starttime=`date +%s`
ip_address=`hostname -I`
ips=($ip_address)
ip=${ips[0]}

thisdir=`pwd`

function test {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "ERROR ${celltype} with $1" >&2
        set_status ${ip} "ERROR" "Error ${celltype} with $1"
        exit 1
    fi
    return $status
}

function set_status {
    local ip=$1
    local stat=$2
    local message=$3
    cd ../progress
    python set-status.py ip:"${ip}" status:"${stat}" message:"${message}"
    cd $thisdir
}

set_status ${ip} "RUNNING" "Started ${celltype} analysis on ${starttime}"

samples=`cat SAMPLES`  # get all arguments
r1_samples=`grep R1_001 SAMPLES`

### Analysis on raw fastq files ###

# FastQC
test ./run-fastqc.sh ${samples}
wait

# Pairwise assembly
set_status ${ip} "RUNNING" "${celltype} Pairwise assembly"
test ./batch-pear.sh ${r1_samples}
wait

### Continue with assembled fastq files ###

samples=`ls *.assembled.fastq.gz`

# Split on MID
set_status ${ip} "RUNNING" "${celltype} Sorting sequences per MID"
test python FastqSplitOnMid.py ${mids} split ${samples}
wait

### Continue with the assembled, split per mid, fastq files ###

samples=`ls split/*.fastq.gz`

# FastQC report
test ./run-fastqc.sh ${samples}
wait

# Search for primers in the fastq files
set_status ${ip} "RUNNING" "${celltype} Searching for primers"
test python motif-search-batch.py ${samples}
wait

# Extract the CDR3 sequence
#set_status ${ip} "RUNNING" "${celltype} Extracting CDR3's"
#test python translate-and-extract-cdr3.py ${celltype} ${samples}
#wait

# Align sequences against IMGT and call SNPs
set_status ${ip} "RUNNING" "${celltype} Aligning sequences"
for ref in $refs; do
    test ./batch-align.sh ${ref} ${samples}
done
wait

### Continue with the aligned sequences ###

bamfiles=`ls *clean.sam`

# Alignment quality report TO IMPLEMENT

## HLAforest for HLA samples
for ref in $refs; do
    for bam in $bamfiles; do
        ./hla-forest.sh ${ref} ${bam}
        wait
    done
done

### Generate reports ###

# Make output directories
mkdir ${beehub_mount}/results-tbcell
mkdir ${beehub_mount}/results-tbcell/raw
mkdir ${beehub_mount}/results-tbcell/reports
mkdir ${beehub_mount}/results-tbcell/final
mkdir ${beehub_mount}/results-tbcell/hla
wait

# Transfer data to Beehub
set_status ${ip} "RUNNING" "Transferring ${celltype} data to Beehub"
test ./copy-to-beehub-hla.sh ${beehub_web}/results-tbcell/hla/
test ./copy-to-beehub-reports.sh ${beehub_web}/results-tbcell/reports/
test ./copy-to-beehub-raw.sh ${beehub_web}/results-tbcell/raw/
cd split
test ./copy-to-beehub-reports.sh ${beehub_web}/results-tbcell/reports/
test ./copy-to-beehub-raw.sh ${beehub_web}/results-tbcell/raw/
cd ../final
test ./copy-to-beehub-reports.sh ${beehub_web}/results-tbcell/reports/
test ./copy-to-beehub-final.sh ${beehub_web}/results-tbcell/final/
cd ..

wait

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"

# Set status when analysis is finished
set_status ${ip} "FINISHED" "${celltype} finished in ${difftime} seconds"
