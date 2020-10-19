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

Manually removed a few files from SAMPLES-fastq. Differences can be seen in the git log

Created an output directory on ResearchDrive for the extracted CDR3 sequences (variables are in ENV.sh):
mkdir /mnt/immunogenomics/RUNS/$RUN/CDR3_20201015

Ran run-translate-and-extract-cdr3.sh to extract the CDR3's and upload the results to the ResearchDrive

Created a notebook to find CDR3s in the files that were produced in the step above: GrepCdr3.ipynb.
Ran this notebook as a test on 2 -CDR3.csv files.
The notebook was exported as: GrepCdr3.py and executed on all -CDR3.csv files.
Input CDR3s as received from Aram: 2020-08-21.presyn-top100-clones-contamination-table-on-pt.csv

The GrepCdr3.py script couldn't read these files:
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/HU5IXOJ02.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run102-r2-HA06QXN02_S2_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run114-r1-HHEXEGX01_S1_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run114-r2-HHEXEGX02_S2_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run114-r4-HHEXEGX04_S4_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run162-r1-HV4JT6A01_S1_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run24-r1-FVUGYZ001_S1_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run37-r3-GAIEXSY03_S3_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run37-r4-GAIEXSY04_S4_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run41-r3-GB5HR3003_S3_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run41-r4-GB5HR3004_S4_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv
ERROR: couldn't read /mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/run58-r2-GMUTUER02_S2_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv

These files are all empty, so that means that no BCRh CDR3's were found in these runs/regions. I'll continue with the remaining runs/regions.

Uploaded the results to: https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/runXXX-roche-LN2/grep-cdr3-20201019/

