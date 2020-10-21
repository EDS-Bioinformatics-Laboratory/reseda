Location of analysis results will end up here (see ENV.sh):
RUN=runXXX-20201021-basta
OUTDIR=autoreactive20201021
https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/runXXX-20201021-basta

Request from Rogier Thurlings (Radboud UMC) to look up sequences that react to autoreactive peptides. Rogier is a collaborator of Guilia Balzarretti.
See e-mails on 9 and 16 October 2020.

Samples:
B005, B007, B008, B009, B010 and B011

Peptides: see file peptides.csv

## Look up in which runs the BASTA samples are

First retrieved all datasheets (metadata) from the researchdrive

```
mkdir Datasheets
cd Datasheets/
ls /mnt/immunogenomics/RUNS/run*miseq/*Datasheet*.csv > DATASHEETS
cp DATASHEETS SAMPLES
../copy-from-beehub.sh 
```

Lookup samples
```
grep -i b005 20* > Basta_Datasheet_20201021.csv
grep -i b007 20* >> Basta_Datasheet_20201021.csv
grep -i b008 20* >> Basta_Datasheet_20201021.csv
grep -i b009 20* >> Basta_Datasheet_20201021.csv
grep -i b010 20* >> Basta_Datasheet_20201021.csv
grep -i b011 20* >> Basta_Datasheet_20201021.csv
```

Manually removed double entries from Basta_Datasheet_20201021.csv (differences can be seen with git diff between two commits)
