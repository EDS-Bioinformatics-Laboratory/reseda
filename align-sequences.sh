#!/bin/bash

ref=$1    # hla_nuc.fasta
#fasta=$2  # sequences.fna
#qual=$3   # sequences.qual
fastq=$2  # sequences.fastq
bwa_param=$3  # optional, e.g.: "-T 20 -B 0"   (don't forget the quotes)

#prefix=`basename ${fasta} .fasta`
mydir=`dirname ${fastq}`
prefix=`basename ${fastq} .fastq.gz`
refprefix=`basename ${ref} .fasta`

#echo "### convert fasta+qual to fastq KEEP FASTQ FILE ###"
#./fasta2fastq.pl -f ${fasta} -q ${qual} -o ${prefix}.fastq
#./fasta2fastq.pl -f ${fasta} -o ${prefix}.fastq
#wait
#gzip -f ${prefix}.fastq
#wait

# tmp directory
mkdir tmp

echo "### align sequences with bwasw ###"
./bwa-0.7.12/bwa mem ${ref} ${mydir}/${prefix}.fastq.gz ${bwa_param} > ${prefix}-${refprefix}.sam
# ./bwa-0.7.12/bwa mem -B 1 -T 20 ${ref} ${mydir}/${prefix}.fastq.gz > ${prefix}-${refprefix}.sam  # keep alignments with lower score
wait

echo "### replace nucleotides that are identical with = ###"
samtools calmd -eS ${prefix}-${refprefix}.sam ${ref} > ${prefix}-${refprefix}-e.sam
wait
rm -f ${prefix}-${refprefix}.sam # REMOVE TMP FILE

echo "### fix CIGAR string KEEP THIS FILE ###"
java -Djava.io.tmpdir=./tmp -jar picard-tools-1.126/picard.jar CleanSam I=${prefix}-${refprefix}-e.sam O=${prefix}-${refprefix}-e-clean.sam
wait
rm -f ${prefix}-${refprefix}-e.sam

echo "### convert sam to bam ###"
java -Djava.io.tmpdir=./tmp -jar ./picard-tools-1.126/picard.jar SamFormatConverter I=${prefix}-${refprefix}-e-clean.sam O=${prefix}-${refprefix}.bam
wait

echo "### add read groups ###"
java -Djava.io.tmpdir=./tmp -jar ./picard-tools-1.126/picard.jar AddOrReplaceReadGroups I=${prefix}-${refprefix}.bam O=${prefix}-${refprefix}RG.bam ID=a LB=a PL=454 PU=MID SM=a SO=coordinate
wait
rm -f ${prefix}-${refprefix}.bam

echo "### index bam file ###"
java -Djava.io.tmpdir=./tmp -jar ./picard-tools-1.126/picard.jar BuildBamIndex I=${prefix}-${refprefix}RG.bam
wait

echo "### Create SAM output that is easy to import in R ###"
samtools view -F 0x04 ${prefix}-${refprefix}RG.bam | cut -f 1,2,3 > ${prefix}-${refprefix}-easy-import.txt

echo "### call raw variants ###"
samtools mpileup -f ${ref} ${prefix}-${refprefix}RG.bam > ${prefix}-${refprefix}RG.pileup
wait
rm -f ${prefix}-${refprefix}RG.bam ${prefix}-${refprefix}RG.bai

if [[ -s ${prefix}-${refprefix}RG.pileup ]]; then
    echo "### call variants that are covered by 1 read ###"
    java -Djava.io.tmpdir=./tmp -jar VarScan.v2.3.7.jar pileup2snp ${prefix}-${refprefix}RG.pileup --min-coverage 1 --min-reads2 1 > ${prefix}-${refprefix}RG.snp.csv
    wait
else
    echo "Pilup file is empty. Skipped VarScan"
fi

rm -f ${prefix}-${refprefix}RG.pileup
wait

echo "### Get unmapped reads ###"
java -Djava.io.tmpdir=./tmp -jar picard-tools-1.126/picard.jar ViewSam I=${prefix}-${refprefix}-e-clean.sam ALIGNMENT_STATUS=Unaligned > ${prefix}-${refprefix}-unmapped.sam
wait

echo "### Unmapped reads sam > fastq ###"
java -Djava.io.tmpdir=./tmp -jar picard-tools-1.126/picard.jar SamToFastq I=${prefix}-${refprefix}-unmapped.sam F=${prefix}-${refprefix}-unmapped.fastq
wait

echo "### Align unmapped reads with less stringent parameters ###"
./bwa-0.7.12/bwa mem ${ref} ${prefix}-${refprefix}-unmapped.fastq -B 0 > ${prefix}-${refprefix}-unmapped-bwa-loose-param.sam
wait

echo "### Only keep the aligned sequences from last alignment ###"
java -Djava.io.tmpdir=./tmp -jar picard-tools-1.126/picard.jar ViewSam I=${prefix}-${refprefix}-unmapped-bwa-loose-param.sam ALIGNMENT_STATUS=Aligned RECORDS_ONLY=true > ${prefix}-${refprefix}-unmapped-aligned-again.sam
wait

rm -f ${prefix}-${refprefix}-unmapped.sam ${prefix}-${refprefix}-unmapped.fastq ${prefix}-${refprefix}-unmapped-bwa-loose-param.sam
