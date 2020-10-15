## Goal

I have received a list of CDR3 sequences from Aram by mail that were shared among several LN samples. We would like to know if these clones come from another sample that is not related to the presyn study.

Approach:

* Run TranslateAndExtractCdr3.py on all raw Roche sequences (input is in fastq format)
* Search for the CDR3 sequences in all runs

Before I can do this I need to convert some raw data from sff to fastq format. In the notebook NOTEBOOKS/ListRawDataRoche.ipynb I compiled lists of runs where the sff files needs to be converted before doing the analysis and I made a list of runs where I already have fastq files.

### Converting sff files to fastq format

I have 4 lists with sff files: SAMPLES-sff-aa SAMPLES-sff-ab SAMPLES-sff-ac SAMPLES-sff-ad

For every list I did the following:

```
cp SAMPLES-sff-aa SAMPLES        # this for every of the 4 files above
nohup ../copy-from-beehub.sh > nohup-get-sff.out 2> nohup-get-sff.err < /dev/null &
bash rename_files.sh
nohup python2 ../SeqToFastq.py sff *.sff > nohup-sfftofastq.out 2> nohup-sfftofastq.err < /dev/null &
```

Uploaded the fastq files with cp-fastq-to-researchdrive-*.sh

### Extract CDR3 sequences from all Roche runs

Listed all the fastq files with NOTEBOOKS/ListRawDataRocheFastq.ipynb

The file SAMPLES-fastq has been created with this notebook
