Data analyzed with git branch runXXX-20201005-maria-reseda20201006

Maria uploaded the zip file HN00134734.zip
I downloaded and unpacked this and renamed the files:

* RenameFilesMaria4.ipynb generated the bash script rename-files.sh
* rename-files.sh really renames files
* The renamed files were uploaded to the "data" directory (for location see ENV.sh)

Created a Datasheet with notebook MakeDataSheet.ipynb

According to Maria the samples 80 and 81 do not have to be analyzed, so I removed it from the Datasheet

Fastq files for the analysis are listed in the file SAMPLES

Ran the analysis in the standalone mode, the variables come from ENV.sh:
mv reference/IGH*human* .
mv mids/MIDS-abc.txt .
mv reftables/ref.table.heavy.csv .
nohup ./execute-all.sh -r $RUN -m MIDS-abc.txt -o $OUTDIR -b no -u no > nohup-exe.out 2> nohup-exe.err < /dev/null &

Made reports and the clones file in the usual way:
python ConcatenateCloneFilesBatch.py -r 20201005-DataSheet-Maria4-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/final/ > tmp.sh
nohup bash tmp.sh > nohup-concat.out 2> nohup-concat.err < /dev/null &
nohup ./report-ALL.sh -r $RUN -i 20201005-DataSheet-Maria4-new.json -b no -o $OUTDIR > nohup-report.out 2> nohup-report.err < /dev/null &
./copy-to-webdav.sh $WEBDAV/$OUTDIR/ cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv assign-info-GC-IGH_HUMAN-after-reassignment.csv vjcdr3-clones-mut-GC-IGH_HUMAN.csv report-*.txt report-all*

Bray-Curtis analysis with SampleSimilarity.ipynb
../copy-to-webdav.sh $WEBDAV/$OUTDIR/similarity/ runMaria4-similarity-IGH_HUMAN*

Look for shared clones with SharedClonesDirection.ipynb
../copy-to-webdav.sh $WEBDAV/$OUTDIR/shared-clones/ runMaria4-IGH-HUMAN-shared-clones*

I can not run RarefactionAnalysis.ipynb in the cloud, so this analysis was done on my laptop

--------------------------

Comparison with the other runs

First renamed the cdr3 file
mv cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv cdr3-clones-GC4-IGH_HUMAN-after-reassignment.csv

Downloaded these cdr3 files from the other runs and renamed them to GC1, GC2, GC3
/mnt/immunogenomics/RUNS/runXXX-20200313-maria/results-tbcell/cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv
/mnt/immunogenomics/RUNS/runXXX-20200420-maria/results-tbcell/cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv
/mnt/immunogenomics/RUNS/runXXX-20200909-maria/results-tbcell/cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv

echo /mnt/immunogenomics/RUNS/runXXX-20200313-maria/results-tbcell/cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv > SAMPLES
../copy-from-beehub.sh 
mv cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv cdr3-clones-GC1-IGH_HUMAN-after-reassignment.csv
echo /mnt/immunogenomics/RUNS/runXXX-20200420-maria/results-tbcell/cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv > SAMPLES
../copy-from-beehub.sh 
mv cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv cdr3-clones-GC2-IGH_HUMAN-after-reassignment.csv
echo /mnt/immunogenomics/RUNS/runXXX-20200909-maria/results-tbcell/cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv > SAMPLES
../copy-from-beehub.sh 
mv cdr3-clones-GC-IGH_HUMAN-after-reassignment.csv cdr3-clones-GC3-IGH_HUMAN-after-reassignment.csv

Got the pt table from the previous run and added the samples of this run (manually)
ls /mnt/immunogenomics/RUNS/runXXX-20200909-maria/20200916-Maria-pt-table-new.csv > SAMPLES
../copy-from-beehub.sh
<manual edit>: 20201007-Maria-pt-table-new.csv
Uploaded this file to the ResearchDrive

Ran the notebook SampleSimilarityAllRuns.ipynb
../copy-to-webdav.sh $WEBDAV/$OUTDIR/similarity/ runMaria1234-similarity-*
