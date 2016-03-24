#!/bin/bash

ref=$1    # hla_nuc.fasta
# fastq=$2   # blah.sff.fastq.gz
# prefix=`basename ${fastq} .fastq.gz`
samfile=$2   # blah.sam
prefix=`basename ${samfile} .sam`

TMP_DIR=`pwd`
mkdir FOREST

# convert fasta+qual to fastq
#./fasta2fastq.pl -f ${fasta} -q ${qual} -o ${prefix}.fastq
#wait
#gzip ${prefix}.fastq
#wait

# align sequences with bwasw
# bwa bwasw -f ${prefix}.sam ${ref} ${prefix}.fastq.gz
# wait

# sort by read name
java -jar picard-tools-1.126/picard.jar SortSam I=${prefix}.sam O=${prefix}.readsort.sam SO=queryname
wait

# Select sequences that are long enough
perl -ne 'if (m/^\@/) {print;} else { @c=split(/\s+/); print if length($c[8]) > 299; }' ${prefix}.readsort.sam > ${prefix}.sizeselection.sam
wait

# in case of multiple hits, choose the first
perl -ne 'if (m/^\@/) {print;} else { @c=split(/\s+/); print if $c[0] ne $previous; $previous=$c[0]; }' ${prefix}.sizeselection.sam > ${prefix}.filtered.sam
wait

# Merge sam files
#java -jar picard-tools-1.126/picard.jar MergeSamFiles I=I4IMECA01.sff.filtered.sam I=I4IMECA02.sff.filtered.sam O=run246.filtered.sam SO=coordinate

# Count alignment hits
perl -ne 'next if m/^@/; @c=split(/\t/); print $c[2],"\n";' ${prefix}.filtered.sam|sort|uniq -c|sort -nr > ${prefix}.filtered.hla.count.txt

# run hlaforest
source config.sh
wait
./build-forest.pl -reads ${prefix}.filtered.sam -o ${prefix} -b $TREES_IN_FOREST $MAX_TREES_FLAG
wait

# scripts below do not work
./prune2.pl -threshold 0.05 -f `echo $TMP_DIR/*.forest | tr " " ","` -t 2 > $TMP_DIR/t2.txt
./prune2.pl -threshold 0.05 -f `echo $TMP_DIR/*.forestpruned | tr " " ","` -t 3 > $TMP_DIR/t3.txt
./prune2.pl -threshold 0.05 -f `echo $TMP_DIR/*.forestprunedpruned | tr " " ","` -t 4 > $TMP_DIR/t4.txt
./prune2.pl -threshold 0.05 -f `echo $TMP_DIR/*.forestprunedprunedpruned | tr " " ","` -t 5 > $TMP_DIR/t5.txt
wait

# call the haplotypes
./call-haplotypes.pl -threshold .3 -2 $TMP_DIR/t2.txt -3 $TMP_DIR/t3.txt -4 $TMP_DIR/t4.txt -5 $TMP_DIR/t5.txt > ${prefix}.haplotypes.txt
wait

# move forest files out of the way
mv $TMP_DIR/*.forest* FOREST
