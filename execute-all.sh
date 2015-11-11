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

# Then run ./execute-all.sh


#######################

starttime=`date +%s`

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

# Split on MID
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

mkdir final

# For each sample; do
for sample in ${samples}; do
    mydir=`dirname ${sample}`
    prefix=`basename ${sample} .fastq.gz`

    # Combine MID, CDR3, V, J and sequence information
    midFile=`echo ${mydir}/${prefix}|perl -ne 's/-.+$/-report.txt/;print;'`
    cdr3File=${sample}-${celltype}-CDR3.csv
    vFile=${prefix}-${v}-easy-import.txt
    jFile=${prefix}-${j}-easy-import.txt
    seqFile=${sample}-${celltype}.csv
    outFile="final/${prefix}-${celltype}-all_info.csv"
    python combine-immuno-data.py ${midFile} ${cdr3File} ${vFile} ${jFile} ${seqFile} ${outFile}
    wait
done

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"
