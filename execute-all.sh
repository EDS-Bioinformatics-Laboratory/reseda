#!/bin/bash

# Preparation:
# mv reference/* .
# mv hlaforest/scripts/* .

# Create a file with all the fastq files:
# ls TESTDATA/* > SAMPLES

# Configure this:
cell="IGH"
organism="human"
celltype="${cell}_HUMAN"
run="runNN-2016MMDD-miseq"

# Reference sequences
mids="MIDS-miseq.txt"
refs="${cell}V_${organism}.fasta ${cell}J_${organism}.fasta"
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
        exit
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
test python fastq-split-on-mid.py ${mids} split ${samples}
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
set_status ${ip} "RUNNING" "${celltype} Extracting CDR3's"
test python translate-and-extract-cdr3.py ${celltype} ${samples}
wait

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
# for ref in $refs; do
#     for bam in $bamfiles; do
#         ./hla-forest.sh ${ref} ${bam}
#         wait
#     done
# done

### Generate reports ###

mkdir final

# For each sample; do
set_status ${ip} "RUNNING" "${celltype} Combining results"
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
    cloneSubsFile="final/${prefix}-${celltype}-clones-subs.csv"
    cloneMainsFile="final/${prefix}-${celltype}-clones-mains.csv"
    totalFile="final/${prefix}-${celltype}-productive.txt"
    test python combine-immuno-data.py ${midFile} ${cdr3File} ${vFile} ${jFile} ${seqFile} ${outFile} ${cloneFile} ${cloneSubsFile} ${cloneMainsFile} ${totalFile}
    wait

done

# Count lines of all_info.csv files
set_status ${ip} "RUNNING" "${celltype} Select correct MIDs"
test wc -l final/*all_info.csv > wc-${ip}.txt
wait
test python select-correct-mids.py wc-${ip}.txt > mv-samples-with-correct-mid.sh
wait
mkdir final/correct-mid
wait
cd final
test bash ../mv-samples-with-correct-mid.sh
wait
mv correct-mid/*-productive.txt .
cd ..

# Correct V gene assignments
set_status ${ip} "RUNNING" "${celltype} Re-assign V genes"
test python re-assign-v-genes.py final/correct-mid/*-all_info.csv
wait
# Move files to final
mv *.rr.* final/correct-mid
wait

# Concatenate clones files
test python concatenate-clone-files.py final/correct-mid/*.rr.clones_subs.csv
wait
mv run-clones_subs.csv run-clones_subs-${ip}.csv
wait

# Make output directories
mkdir ${beehub_mount}/results-tbcell
mkdir ${beehub_mount}/results-tbcell/raw
mkdir ${beehub_mount}/results-tbcell/reports
mkdir ${beehub_mount}/results-tbcell/final
mkdir ${beehub_mount}/results-tbcell/final/correct-mid
mkdir ${beehub_mount}/results-tbcell/hla
wait

# Transfer data to Beehub
set_status ${ip} "RUNNING" "Transferring ${celltype} data to Beehub"
test curl -T run-clones_subs-${ip}.csv --netrc ${beehub_web}/results-tbcell/
test ./copy-to-beehub-reports.sh ${beehub_web}/results-tbcell/reports/
test ./copy-to-beehub-raw.sh ${beehub_web}/results-tbcell/raw/
test ./copy-to-beehub-hla.sh ${beehub_web}/results-tbcell/hla/
cd split
test ./copy-to-beehub-reports.sh ${beehub_web}/results-tbcell/reports/
test ./copy-to-beehub-raw.sh ${beehub_web}/results-tbcell/raw/
cd ../final
test ./copy-to-beehub-reports.sh ${beehub_web}/results-tbcell/reports/
test ./copy-to-beehub-final.sh ${beehub_web}/results-tbcell/final/
cd correct-mid
test ./copy-to-beehub-final.sh ${beehub_web}/results-tbcell/final/correct-mid/
cd ../..

wait

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"

# Set status when analysis is finished
set_status ${ip} "FINISHED" "${celltype} finished in ${difftime} seconds"
