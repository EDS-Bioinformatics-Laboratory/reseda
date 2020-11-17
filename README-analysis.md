## Data

RUN=run40-20201116-miseq
OUTDIR=reseda20201116

Data location: https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/run40-20201116-miseq

## Data analysis

Data analyzed with git branch run40-20201116-miseq-reseda20201116

The RESEDA analysis was done according to the instructions in the documentation / README.md file: "How to run - Using PiCaS for sending jobs and using a job monitor tool"

```
878  source ENV.sh 
  879  ls /mnt/immunogenomics/RUNS/$RUN/metadata/* > SAMPLES 
  880  cat SAMPLES 
  881  nano SAMPLES 
  882  ./copy-from-beehub.sh 
  883  python VerifyBasespaceCopy.py 
  884  python MetaData.py 20201116-RUN40-DataSheet.csv > 20201116-RUN40-DataSheet.json
  885  python VerifyBasespaceCopy.py -i 20201116-RUN40-DataSheet.json -r $RUN 
  886  nohup bash basespace-calc-checksum.sh > nohup-checksum.out 2> nohup-checksum.err < /dev/null &
  887  lr
  888  python MetaData.py 20201116-RUN40-DataSheet-no-umi.csv > 20201116-RUN40-DataSheet-no-umi.json
  889  python MetaData.py 20201116-RUN40-DataSheet-UMI-Constant.csv > 20201116-RUN40-DataSheet-UMI-Constant.json
  890  python MetaData.py 20201116-RUN40-DataSheet-UMI-Constant-RACE.csv > 20201116-RUN40-DataSheet-UMI-Constant-RACE.json
  891  python MetaData.py 20201116-RUN40-DataSheet-UMI.csv > 20201116-RUN40-DataSheet-UMI.json
  892  python VerifyBasespaceCopy.py -i 20201116-RUN40-DataSheet.json -r $RUN 
  893  less VerifyBasespaceCopy.py 
  894  lr
  895  ./copy-to-webdav.sh $WEBDAV/metadata/ *.json
  896  ls
  897  mkdir tokens
  898  python 
  899  python MakeSamplesFiles.py 
  900  ls *.json
  901  python MakeSamplesFiles.py -r 20201116-RUN40-DataSheet-no-umi.json -w /mnt/immunogenomics/RUNS/$RUN/data/
  902  lr
  903  wc -l SAMPLES-run40-mouse-BCRh
  904  mv SAMPLES-run40-mouse-BCRh SAMPLES-run40-no-umi-mouse-BCRh
  905  python MakeSamplesFiles.py -r 20201116-RUN40-DataSheet-UMI.json -w /mnt/immunogenomics/RUNS/$RUN/data/
  906  lr
  907  mv SAMPLES-run40-human-TCRb SAMPLES-run40-UMI-human-TCRb
  908  mv SAMPLES-run40-human-BCRh SAMPLES-run40-UMI-human-BCRh
  909  python MakeSamplesFiles.py -r 20201116-RUN40-DataSheet-UMI-Constant.json -w /mnt/immunogenomics/RUNS/$RUN/data/
  910  lr
  911  mv SAMPLES-run40-human-BCRh SAMPLES-run40-UMI-Constant-human-BCRh
  912  python MakeSamplesFiles.py -r 20201116-RUN40-DataSheet-UMI-Constant-RACE.json -w /mnt/immunogenomics/RUNS/$RUN/data/
  913  lr
  914  mv SAMPLES-run40-human-BCRh SAMPLES-run40-UMI-Constant-RACE-human-BCRh
```

## After the analysis

Done according the documentation. The RACE samples are not analyzed yet (17 Nov 2020), because a few changes are needed in the scripts and in the MID definition file

* ConcatenateCloneFilesBatch.py
* report-ALL.sh
* SampleSimilarity (notebook in repository)
* SharedClonesDirection.ipynb (notebook in repository)
