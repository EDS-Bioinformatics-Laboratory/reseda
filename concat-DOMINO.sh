RUN=project-20210304-DOMINO
OUTDIR=processing/20210304-reseda

python ConcatenateCloneFilesBatch.py -r 20170821-RUN17-datasheet-UMI-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run17/final/ > concat.sh
python ConcatenateCloneFilesBatch.py -r 20171005_RUN019_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run19/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20171127_RUN20-Datasheet-UMI-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run20/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20171218_RUN021_Datasheet-umi-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run21/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20180119_RUN22_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run22/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20180302_RUN23_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run23/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20180409_RUN26_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run26/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20180514_RUN28_Datasheet-UMI-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run28/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20180605_RUN29_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run29/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20180723_RUN30_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run30/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20180921_RUN31_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run31/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20181125_RUN32_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run32/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20181214_RUN33_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run33/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20190214_RUN34_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run34/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20190609_RUN35_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run35/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20191003_RUN36_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run36/final/ >> concat.sh
python ConcatenateCloneFilesBatch.py -r 20200529_RUN38_Datasheet-DOMINO-new.json -w /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run38/final/ >> concat.sh

