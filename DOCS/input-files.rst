Prepare input files
===================

Reference database
------------------

The reference files are already part of the repository. The following steps are necessary for an update of the reference database or for a new reference.

**Required for V and J assignment**

Download the reference sequences of the chains you are interested in (nucleotide sequences in fasta format) from the IMGT website (http://imgt.org/)

* Select species, gene type, functionality (functional) and the locus in the GeneDB
* Select all genes and download the "F+ORF+all P nucleotide sequences", store as e.g. TRBV_human.fasta and TRBJ_human.fasta
* Do this for the variable and joining genes
* Build a BWA index and a samtools faidx on the fasta files

**Required for CDR3 identification**

Download the peptide sequences of the variable genes (F+ORF+in-frame P amino acid sequences with IMGT gaps). Convert the downloaded fasta entries with the helper-ref-table.py script

**Adding reference sequences of a new organism**

The new reference files need to be specified in TranslateAndExtractCdr3.py

Expected input (sequences)
--------------------------

RESEDA assumes that input files have the following names and comes in fastq format:

Paired end: SAMPLENAME_S1_L001_R1_001.fastq.gz and SAMPLENAME_S1_L001_R2_001.fastq.gz

Where SAMPLENAME is unique in a run and S1 can have any number.

Single end: SAMPLENAME_S1_L001.assembled.fastq.gz

Data descriptions (metadata)
----------------------------

When the data analysis is done with ToPoS you need a 'Datasheet' from which input parameters are generated automatically.
An example can be found in the 'TESTDATA' directory: 20190609_RUN35_Datasheet-UMI.csv

The columns 'Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2' can be left empty.
If there are both sample with and without UMI in the sequence run the datasheet needs to be split into one datasheet with all the samples with and one datasheet with the samples without UMI.
It is easiest to create different output directories for these two types of samples, but this is not mandatory.

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
