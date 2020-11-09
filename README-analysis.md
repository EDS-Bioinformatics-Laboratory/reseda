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

## Reports

Made a datasheet with ``MakeDataSheet.ipynb``: 20201013-DataSheet-Maria-mouse.csv

Converted the datasheet to json format as I always do and uploaded the metadata to the researchdrive

```
python MetaData.py 20201013-DataSheet-Maria-mouse.csv > 20201013-DataSheet-Maria-mouse.json
python MakeSamplesFiles.py -r 20201013-DataSheet-Maria-mouse.json -w /mnt/immunogenomics/RUNS/$RUN/data/
./copy-to-webdav.sh $WEBDAV/ 20201013-DataSheet-Maria-mouse*
```

Concatenated the clones files and made the run report, the usual way

```
python ConcatenateCloneFilesBatch.py -r 20201013-DataSheet-Maria-mouse-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/final/ > tmp.sh
nohup bash tmp.sh > nohup-concat.out 2> nohup-concat.err < /dev/null &
nohup ./report-ALL.sh -r $RUN -i 20201013-DataSheet-Maria-mouse-new.json -b no -o $OUTDIR > nohup-report.out 2> nohup-report.err < /dev/null &
./copy-to-webdav.sh $WEBDAV/$OUTDIR/ report-*.txt report-all* cdr3-clones-GC-IGH_MOUSE-after-reassignment.csv assign-info-GC-IGH_MOUSE-after-reassignment.csv vjcdr3-clones-mut-GC-IGH_MOUSE.csv
```

## Notebooks

Moved files to the NOTEBOOKS directory

```
mv 20201013-DataSheet-Maria-mouse.csv NOTEBOOKS/
mv cdr3-clones-GC-IGH_MOUSE-after-reassignment.csv NOTEBOOKS/
mv vjcdr3-clones-mut-GC-IGH_MOUSE.csv NOTEBOOKS/
```

Ran:

* SampleSimilarity.ipynb
* SharedClonesDirection.ipynb
