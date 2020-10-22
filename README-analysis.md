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

Manually removed the Behcet samples from the list.

I couldn't find all samples, so I asked Rogier for the alternative names of these samples. Answer:

B005      0352      SS3
B007      4708
B008      2172
B009      1193      SS4
B010      3684
B011      8296      1467     2018mrt383

Made a new list with the samples based on the information above:

```
grep -i basta 20*|grep -i b005 > Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i 0352 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i SS3 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i B007 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i 4708 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i B008 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i 2172 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i B009 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i 1193 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i SS4 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i B010 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i 3684 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i B011 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i 8296 >> Basta_Datasheet_20201022.csv 
grep -i basta 20*|grep -i 1467 >> Basta_Datasheet_20201022.csv 
```

Made the list uniq

``perl -ne '@c=split(/:/); print $c[$#c];' Basta_Datasheet_20201022.csv |sort|uniq > Basta_Datasheet_20201022_uniq.csv``

