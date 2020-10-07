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

