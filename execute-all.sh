#!/bin/bash

# Preparation:
# mv reference/* .
# mv hlaforest/scripts/* .

# Create a file with all the fastq files
# ls TESTDATA/* > SAMPLES

starttime=`date +%s`

celltype="IGH_HUMAN"
mids="MIDS-miseq.txt"
refs="IGHV_human.fasta IGHJ_human.fasta"
samples=`cat SAMPLES`  # get all arguments

r1_samples=`grep R1_001 SAMPLES`

### Analysis on raw fastq files ###

# FastQC
./run-fastqc.sh ${samples}
wait

# Pairwise assembly
./batch-pear.sh ${r1_samples}
wait

### Continue with assembled fastq files ###

samples=`ls *.assembled.fastq.gz`

# Split on MID (TO MODIFY)
python fastq-split-on-mid.py ${mids} split ${samples}
wait

### Continue with the assembled, split per mid, fastq files ###

samples=`ls split/*.fastq.gz`

# FastQC report
./run-fastqc.sh ${samples}
wait

# Search for primers in the fastq files
python motif-search-batch.py ${samples}
wait

# Extract the CDR3 sequence
python translate-and-extract-cdr3.py ${celltype} ${samples}
wait

# Align sequences against IMGT and call SNPs
for ref in $refs; do
    ./batch-align.sh ${ref} ${samples}
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

# For each sample; do

# Combine MID, CDR3, V, J and sequence information
# midFile="/mnt/immunogenomics/RUNS/run03-20150814-miseq/data/midsort/BCRh_S40_L001_R2_001-report.txt"
# cdr3File="/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/paul/BCRh_S40_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv"
# vFile="/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/BCRh/BCRh_S40_L001.assembled-IGHV_human-easy-import.txt"
# jFile="/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/BCRh/BCRh_S40_L001.assembled-IGHJ_human-easy-import.txt"
# seqFile="/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/paul/BCRh_S40_L001.assembled.fastq.gz-IGH_HUMAN.csv"
# outFile="all_info.txt"
# python combine-immuno-data.py ${midFile} ${cdr3File} ${vFile} ${jFile} ${seqFile} ${outFile}
# wait

# done

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"
