Prepare input files
===================

Expected input
--------------

RESEDA assumes that input files have the following names and comes in fastq format:

Paired end: SAMPLENAME_S1_L001_R1_001.fastq.gz and SAMPLENAME_S1_L001_R2_001.fastq.gz
Where SAMPLENAME is unique in a run and S1 can have any number.

Single end: SAMPLENAME_S1_L001.assembled.fastq.gz

Prepare Roche data
------------------

* Use MakePTtableFromAAreads.R - Create a pt.table (sample description) from AA.reads file
* SplitAAreads.py - Splits the AA.reads table per sample (check the column names that you want to include in the file name!)
* SeqToFastq.py - Convert the tab-delimited files to fastq format

When starting from files in sff format:

* Convert sff to fastq with SeqToFastq.py
* Split the fastq files with FastqSplitOnMid.py
* Rename the files with jupyter notebook RenameSplitRoche.ipynb
