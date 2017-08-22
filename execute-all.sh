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
protocol=$6

# Or configure settings here:
# run="runNN-00-08-20160722-human-BCRh"
# mids="MIDS-miseq.txt"
# organism="human"
# cell="IGH"
# celltype="${cell}_HUMAN"
# protocol="paired"          # "paired" or "single" end (default: paired)

# Reference sequences
refs="${cell}V_${organism}.fasta ${cell}J_${organism}.fasta"
v="${cell}V_${organism}"
j="${cell}J_${organism}"

# Mount the Beehub webdav server and configure the location
resultsdir="results-tbcell"
beehub_mount="/mnt/immunogenomics/RUNS/${run}"
beehub_web="https://beehub.nl/amc-immunogenomics/RUNS/${run}"

# Then run ./execute-all.sh


#######################

starttime=`date +%s`
ip_address=`hostname -I`
ips=($ip_address)
ip=${ips[0]}

thisdir=`pwd`

# Paired-end sequences are the default
if [ ${protocol} -ne "single" ]; then
    protocol="paired"
fi

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

# Job monitoring. To use this you need https://bitbucket.org/barbera/progress
function set_status {
    local ip=$1
    local stat=$2
    local message=$3
    cd ../progress
    python set-status.py ip:"${ip}" status:"${stat}" message:"${message}"
    cd $thisdir
}

runcmd ./log-versions.sh

set_status ${ip} "RUNNING" "Started ${mids} ${celltype} analysis on ${starttime}"

# Remote files:
runcmd ./copy-from-beehub.sh
wait
# # Local files:
# cp SAMPLES LOCAL_SAMPLES
# wait

samples=`cat LOCAL_SAMPLES`  # get all arguments
r1_samples=`grep R1_001 LOCAL_SAMPLES`

### Analysis on raw fastq files ###
if [[ ${protocol} -eq "paired" ]]; then
    # FastQC
    runcmd ./run-fastqc.sh ${samples}
    wait

    # Pairwise assembly
    set_status ${ip} "RUNNING" "${celltype} Pairwise assembly"
    runcmd ./batch-pear.sh ${r1_samples}
    wait

    samples=`ls *.assembled.fastq.gz`
fi

### Continue with assembled fastq files ###

# Split on MID
set_status ${ip} "RUNNING" "${celltype} Sorting sequences per MID"
runcmd python2 FastqSplitOnMid.py ${mids} split ${samples}
wait

### Continue with the assembled, split per mid, fastq files ###

samples=`ls split/*.fastq.gz`

# FastQC report
runcmd ./run-fastqc.sh ${samples}
wait

# # Search for primers in the fastq files
# set_status ${ip} "RUNNING" "${celltype} Searching for primers"
# runcmd python2 motif-search-batch.py ${samples}
# wait

# Extract the CDR3 sequence
set_status ${ip} "RUNNING" "${celltype} Extracting CDR3's"
runcmd python2 TranslateAndExtractCdr3.py -c ${celltype} ${samples}
wait

# Align sequences against IMGT and call SNPs
set_status ${ip} "RUNNING" "${celltype} Aligning sequences"
for ref in $refs; do
    runcmd ./batch-align.sh ${ref} ${samples}
done
wait

set_status ${ip} "RUNNING" "Determine SNPs from the SAM files" # creates file: ${prefix}-${refprefix}-e-clean.sam.mut.txt
runcmd python MutationsFromSam.py *-e-clean.sam
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
    runcmd python2 combine-immuno-data.py ${midFile} ${cdr3File} ${vFile} ${jFile} ${seqFile} ${outFile} ${cloneFile} ${cloneSubsFile} ${cloneMainsFile} ${totalFile}
    wait

done

# Count lines of all_info.csv files
set_status ${ip} "RUNNING" "${celltype} Select correct MIDs"
runcmd wc -l final/*all_info.csv > wc-${ip}.txt
wait
runcmd python2 select-correct-mids.py wc-${ip}.txt > mv-samples-with-correct-mid.sh
wait
mkdir final/correct-mid
wait
cd final
runcmd bash ../mv-samples-with-correct-mid.sh
wait
mv correct-mid/*-productive.txt .
cd ..

# Correct V gene assignments
set_status ${ip} "RUNNING" "${celltype} Re-assign V genes"
runcmd python2 re-assign-v-genes.py final/correct-mid/*-all_info.csv
wait
# Move files to final
mv *.rr.* final/correct-mid
wait

# Do mutation analysis if reference is IGH_HUMAN
if [[ ${celltype} -eq "IGH_HUMAN" ]]; then
    set_status ${ip} "RUNNING" "${celltype} Mutation analysis"
    samples=`ls final/correct-mid/*-IGH_HUMAN-all_info.csv.rr.all_info.csv`
    for sample in ${samples}; do
        prefix=`basename ${sample} -IGH_HUMAN-all_info.csv.rr.all_info.csv`
        runcmd Rscript MutationAnalysisVJ.R indir=\'.\' outdir=\'final/correct-mid\' V.file=\'${prefix}-IGHV_human-e-clean.sam.mut.txt\' J.file=\'${prefix}-IGHJ_human-e-clean.sam.mut.txt\' CDR3.file=\'final/correct-mid/${prefix}-IGH_HUMAN-all_info.csv.rr.all_info.csv\'
    done
fi

# Make output directories
mkdir ${beehub_mount}
mkdir ${beehub_mount}/${resultsdir}
mkdir ${beehub_mount}/${resultsdir}/raw
mkdir ${beehub_mount}/${resultsdir}/raw/correct-mid
mkdir ${beehub_mount}/${resultsdir}/reports
mkdir ${beehub_mount}/${resultsdir}/final
mkdir ${beehub_mount}/${resultsdir}/final/correct-mid
mkdir ${beehub_mount}/${resultsdir}/hla
wait

# Transfer data to Beehub
set_status ${ip} "RUNNING" "Transferring ${celltype} data to Beehub"
#runcmd curl -T run-clones_subs-${ip}.csv --netrc ${beehub_web}/${resultsdir}/
runcmd ./copy-to-beehub-reports.sh ${beehub_web}/${resultsdir}/reports/
runcmd ./copy-to-beehub-raw.sh ${beehub_web}/${resultsdir}/raw/
runcmd ./copy-to-beehub-hla.sh ${beehub_web}/${resultsdir}/hla/
cd split
runcmd ./copy-to-beehub-reports.sh ${beehub_web}/${resultsdir}/reports/
runcmd ./copy-to-beehub-raw.sh ${beehub_web}/${resultsdir}/raw/
cd ../final
runcmd ./copy-to-beehub-reports.sh ${beehub_web}/${resultsdir}/reports/
runcmd ./copy-to-beehub-final.sh ${beehub_web}/${resultsdir}/final/
cd correct-mid
runcmd ./copy-to-beehub-final.sh ${beehub_web}/${resultsdir}/final/correct-mid/
runcmd ./copy-to-beehub-final.sh ${beehub_web}/${resultsdir}/raw/correct-mid/
cd ../..

wait

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"

# Set status when analysis is finished
set_status ${ip} "FINISHED" "${celltype} finished in ${difftime} seconds"
