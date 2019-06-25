Prepare input files
===================

Reference database
------------------

**Required for V and J assignment**

Download the reference sequences of the chains you are interested in (nucleotide sequences in fasta format) from the IMGT website (http://imgt.org/)

* Select species, gene type, functionality (functional) and the locus in the GeneDB
* Select all genes and download the "F+ORF+all P nucleotide sequences", store as e.g. TRBV_human.fasta and TRBJ_human.fasta
* Do this for the variable and joining genes
* Build a BWA index and a samtools faidx on the fasta files

**Required for CDR3 identification**

Download the peptide sequences of the variable genes (F+ORF+in-frame P amino acid sequences with IMGT gaps). Convert the downloaded fasta entries with the helper-ref-table.py script

Expected input
--------------

RESEDA assumes that input files have the following names and comes in fastq format:

Paired end: SAMPLENAME_S1_L001_R1_001.fastq.gz and SAMPLENAME_S1_L001_R2_001.fastq.gz

Where SAMPLENAME is unique in a run and S1 can have any number.

Single end: SAMPLENAME_S1_L001.assembled.fastq.gz

Prepare Roche data
------------------

When starting from a table with sequence data. Note that this table usually
doesn't contain quality scores. These will be set to 40.

* Use MakePTtableFromAAreads.R - Create a pt.table (sample description) from AA.reads file
* SplitAAreads.py - Splits the AA.reads table per sample (check the column names that you want to include in the file name!)
* SeqToFastq.py - Convert the tab-delimited files to fastq format

When starting from files in sff format:

* Convert sff to fastq with SeqToFastq.py
* Split the fastq files with FastqSplitOnMid.py
* Rename the files with jupyter notebook RenameSplitRoche.ipynb
