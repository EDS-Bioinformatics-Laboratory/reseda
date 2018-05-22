#!/bin/bash

# Preparation:
# mv reference/* .
# mv hlaforest/scripts/* .

# Create a file with all the fastq files:
# ls TESTDATA/* > SAMPLES

# Get settings from commandline arguments
#run=$1
#mids=$2
#organism=$3
#cell=$4
#celltype=$5

# Configure this:
run="run28-20180514-miseq"
mids="MIDS-miseq.txt"
cell="HLA"
organism="human"
celltype="${cell}_HUMAN"

# Reference sequences
refs="hla_nuc_nospace.fasta"

# Mount the Beehub webdav server and configure the location
resultsdir="hla"
beehub_mount="/mnt/immunogenomics/RUNS/${run}"
beehub_web="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/${run}"

# Then run ./execute-all-hla.sh


#######################

starttime=`date +%s`
ip_address=`hostname -I`
ips=($ip_address)
ip=${ips[0]}

thisdir=`pwd`

function runcmd {
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
runcmd ./run-fastqc.sh ${samples}
wait

# Pairwise assembly
set_status ${ip} "RUNNING" "${celltype} Pairwise assembly"
runcmd ./batch-pear.sh ${r1_samples}
wait

### Continue with assembled fastq files ###

samples=`ls *.assembled.fastq.gz`

# # Split on sequence length
# set_status ${ip} "RUNNING" "${CELLTYPE} Split sequences on length"
# runcmd python2 FastqSplitOnSequenceLength.py -l 270 ${samples}
# wait
#
# # New sample list
# samples=`cat SAMPLES_long`

# Split on MID
set_status ${ip} "RUNNING" "${celltype} Sorting sequences per MID"
runcmd python FastqSplitOnMid.py no ${mids} split ${samples}
wait

### Continue with the assembled, split per mid, fastq files ###

samples=`ls split/*.fastq.gz`

# FastQC report
runcmd ./run-fastqc.sh ${samples}
wait

# Align sequences against IMGT and call SNPs
set_status ${ip} "RUNNING" "${celltype} Aligning sequences"
for ref in $refs; do
    runcmd ./batch-align.sh ${ref} ${samples}
done
wait

### Continue with the aligned sequences ###

bamfiles=`ls *clean.sam`

# Alignment quality report TO IMPLEMENT

# HLAforest for HLA samples
for ref in $refs; do
    for bam in $bamfiles; do
        ./hla-forest.sh ${ref} ${bam}
        wait
    done
done

### Generate reports ###

# Make output directories
mkdir ${beehub_mount}
mkdir ${beehub_mount}/${resultsdir}
mkdir ${beehub_mount}/${resultsdir}/raw
mkdir ${beehub_mount}/${resultsdir}/reports
mkdir ${beehub_mount}/${resultsdir}/final
mkdir ${beehub_mount}/${resultsdir}/final/correct-mid
mkdir ${beehub_mount}/${resultsdir}/hla
wait

# Transfer data to Beehub
set_status ${ip} "RUNNING" "Transferring ${CELLTYPE} data to Webdav"
runcmd ./copy-to-webdav.sh ${beehub_web}/${resultsdir}/reports/ *-pear.log *-pear.err *.quality-filter.log wc-*.txt versions-*
runcmd ./copy-to-webdav.sh ${beehub_web}/${resultsdir}/raw/ *.sam *.snp.csv *.mut.txt *.short*.assembled.fastq.gz
runcmd ./copy-to-webdav.sh ${beehub_web}/${resultsdir}/hla/ *.hla.count.txt* *.haplotypes.txt *.seqlength.report

runcmd ./copy-to-webdav.sh ${beehub_web}/${resultsdir}/reports/ split/*.primers.count.txt split/*-report.txt split/*-midcount.txt split/*-extra.txt
runcmd ./copy-to-webdav.sh ${beehub_web}/${resultsdir}/raw/ split/*.fastq.gz split/*_fastqc.zip split/*-alt-V-CDR3.csv split/*-alt-J-CDR3.csv

runcmd ./copy-to-webdav.sh ${beehub_web}/${resultsdir}/reports/ final/*-productive.txt
runcmd ./copy-to-webdav.sh ${beehub_web}/${resultsdir}/final/ final/*-all_info.csv final/*-clones-subs.csv

runcmd ./copy-to-webdav.sh ${beehub_web}/${resultsdir}/final/correct-mid/ final/correct-mid/*.rr.* final/correct-mid/*mutations*
runcmd ./copy-to-webdav.sh ${beehub_web}/${resultsdir}/raw/correct-mid/ final/correct-mid/*-all_info.csv final/correct-mid/*-clones-subs.csv

wait

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"

# Set status when analysis is finished
set_status ${ip} "FINISHED" "${celltype} finished in ${difftime} seconds"
