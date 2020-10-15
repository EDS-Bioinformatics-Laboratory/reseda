myfiles=`cat NOTEBOOKS/SAMPLES-fastq`

for myfile in ${myfiles}; do
    # Get the file from the ResearchDrive
    echo ${myfile}
    echo ${myfile} > SAMPLES
    ./copy-from-beehub.sh
    wait

    # Run TranslateAndExtractCdr3.py
    fastqfile=`cat LOCAL_SAMPLES`
    python2 TranslateAndExtractCdr3.py ${fastqfile}
    wait

    # Upload the results and remove the local copy
    ./copy-to-webdav.sh "https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/" ${fastqfile}-IGH_HUMAN-CDR3.csv ${fastqfile}-IGH_HUMAN-report.txt
    wait

    # Remove local copies
    rm ${fastqfile}*
    wait
done

echo "FINISHED"
