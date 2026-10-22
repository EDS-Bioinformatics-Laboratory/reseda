#!/bin/bash

RUN=$1
if [ -z "$RUN" ]; then
    echo "Usage: $0 <RUN>"
    exit 1
fi

prefix="remote"

# rclone ls remote:RUNS/project-20261014-cea-cef/Data/NameOfDataset_1/Raw
remotepath="$prefix:RUNS/${RUN}/Data/NameOfDataset_1/Raw"
echo "Remote path: ${remotepath}"

starttime=`date +%s`

samples=`perl -ne 's/.* //;print;' < SAMPLES`  # get all samples: /mnt/immunogenomics/RUNS/runXXX/Data/NameOfDataset_1/Raw/*.fastq.gz

echo "Transferring files..."
rm -f LOCAL_SAMPLES
for sample in $samples; do
    samplefile=`echo ${remotepath}/${sample}`
    echo ${samplefile}
    localsample=`basename ${sample}`
    rclone copyto ${samplefile} ${localsample} --progress --transfers=4 --checkers=8 --retries=3
    echo ${localsample} >> LOCAL_SAMPLES
    wait
done

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`

echo "Finished in $difftime seconds"

# Transfer files FROM beehub (-C - is for resuming files)
# curl -C - -L -u bschaik -O https://beehub.nl/home/bschaik/immunology/RUNS/run234/data/blah.sff
