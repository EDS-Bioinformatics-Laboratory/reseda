RUN=run61
WEBDAV=https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/${RUN}/data/
../copy-to-webdav.sh $WEBDAV ${RUN}-*.fastq.gz

RUN=run70
WEBDAV=https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/${RUN}/data/
../copy-to-webdav.sh $WEBDAV ${RUN}-*.fastq.gz

RUN=run71
WEBDAV=https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/${RUN}/data/
../copy-to-webdav.sh $WEBDAV ${RUN}-*.fastq.gz

echo "FINISHED"
