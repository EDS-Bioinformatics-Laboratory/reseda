## Data analysis mouse GC data

Collaborators:

* Maria
* Jeroen Guikema

RESEDA analysis:

* Barbera van Schaik
* Antoine van Kampen

## Source

Maria uploaded raw data to https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/runXXX-20201013-maria-mouse/HN00133040/raw-data/

Barbera downloaded this dataset for the RESEDA analysis

Raw data was renamed with ``RenameFastq.ipynb`` and ``rename.sh``.
The result was uploaded to: https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/runXXX-20201013-maria-mouse/data/

```
source ENV.sh   # This file contains location on the researchdrive
bash rename.sh 
./copy-to-webdav.sh $WEBDAV/data/ *.fastq.gz
```

## Analysis

Prepare analysis (moved reference fastq and table to working directory)

```
mv reference/IGH*mouse* .
mv mids/MIDS-abc.txt .
mv reftables/ref.table.mouse.heavy.csv .
```

Started the analysis in the stand-alone mode (see documentation)

```
mkdir maria-mouse
mv *.fastq.gz maria-mouse/
ls maria-mouse/* > SAMPLES 
nohup ./execute-all.sh -r $RUN -l local -m MIDS-abc.txt -org mouse -celltype IGH_MOUSE -o $OUTDIR -b no -u no > nohup-exe.out 2> nohup-exe.err < /dev/null &
```

