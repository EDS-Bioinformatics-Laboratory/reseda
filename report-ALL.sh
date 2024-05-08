#!/bin/bash

################## Get arguments #################

function show_help {
    echo "USAGE: ./report-ALL.sh [options]"
    echo "    -r --run             required: runNN-YYYYMMDD-miseq"
    echo "    -i --info            required: yyyymmdd_RUNnn_Datasheet-new.json"
    echo "    -b --barcodes        default: yes (other option: no)"
    echo "    -o --outdir          default: results-tbcell"
    exit
}

OUTDIR="results-tbcell"
BARCODES="yes"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -h|--help)
        show_help
        shift # past argument
        ;;
        -r|--run)
        RUN="$2"
        shift # past argument
        shift # past value
        ;;
        -o|--outdir)
        OUTDIR="$2"
        shift # past argument
        shift # past value
        ;;
        -i|--info)
        INFO="$2"
        shift # past argument
        shift # past value
        ;;
        -b|--barcodes)
        BARCODES="$2"
        shift # past argument
        shift # past value
        ;;
        *)    # unknown option
        POSITIONAL+=("$1") # save it in an array for later
        shift # past argument
        ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

# Check input
if [ "${RUN}" == "" ]; then
    show_help
fi
if [ "${INFO}" == "" ]; then
    show_help
fi

echo RUN             = "${RUN}"
echo INFO            = "${INFO}"
echo OUTDIR          = "${OUTDIR}"
echo BARCODES        = "${BARCODES}"
# echo POSITIONAL     = "${POSITIONAL}"
# echo rest           = "$@"

################## End get arguments #############

files="/mnt/immunogenomics/RUNS/${RUN}/Processing/${OUTDIR}/Results/reports/*-pear.log"
wait
echo $files > SAMPLES
./copy-from-webdav.sh
wait
localfiles=`cat LOCAL_SAMPLES`
grep '^Assembled reads \.' $localfiles > report-PEAR.txt
wait
rm -f $localfiles
echo "Wrote report-PEAR.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/Processing/${OUTDIR}/Results/reports/*-midcount.txt"
wait
echo $files > SAMPLES
./copy-from-webdav.sh
wait
localfiles=`cat LOCAL_SAMPLES`
more $localfiles > report-MIDs.txt
wait
rm -f $localfiles
echo "Wrote report-MIDs.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/Processing/${OUTDIR}/Results/reports/*HUMAN-report.txt /mnt/immunogenomics/RUNS/${RUN}/Processing/${OUTDIR}/Results/reports/*MOUSE-report.txt /mnt/immunogenomics/RUNS/${RUN}/Processing/${OUTDIR}/Results/reports/*RHESUS-report.txt"
wait
echo $files > SAMPLES
./copy-from-webdav.sh
wait
localfiles=`cat LOCAL_SAMPLES`
grep '^4' $localfiles > report-CDR3.txt
echo "Wrote report-CDR3.txt"
grep '^5' $localfiles > report-ALT-V.txt
echo "Wrote report-ALT-V.txt"
grep '^6' $localfiles > report-ALT-J.txt
wait
rm -f $localfiles
echo "Wrote report-ALT-J.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/Processing/${OUTDIR}/Results/reports/*-productive.txt"
wait
echo $files > SAMPLES
./copy-from-webdav.sh
wait
localfiles=`cat LOCAL_SAMPLES`
grep 'Unique reads with V and J in all_info' $localfiles > report-PRODUCTIVE.txt
wait
rm -f $localfiles
echo "Wrote report-PRODUCTIVE.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/Processing/${OUTDIR}/Results/raw/*.sam"
wait
echo $files > SAMPLES
./copy-from-webdav.sh
wait
localfiles=`cat LOCAL_SAMPLES`
python report-alignments.py ${localfiles}
wait
rm -f $localfiles
echo "Wrote report-ALIGNED-*"

files="/mnt/immunogenomics/RUNS/${RUN}/Processing/${OUTDIR}/Results/final/*.rr.clones_subs.csv"
wait
echo $files > SAMPLES
./copy-from-webdav.sh
wait
localfiles=`cat LOCAL_SAMPLES`
python2 report-after-v-reassignment.py $localfiles
wait
rm -f $localfiles
echo "Wrote report-AFTER-V-REASSIGNMENT.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/Processing/${OUTDIR}/Results/reports/*.quality-filter.log"
wait
echo $files > SAMPLES
./copy-from-webdav.sh
wait
localfiles=`cat LOCAL_SAMPLES`
perl -ne 'print if ! m/^datafile/;' $localfiles > report-QUALITY-FILTER.txt
wait
rm -f $localfiles
echo "Wrote report-QUALITY-FILTER.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/Processing/${OUTDIR}/Results/reports/*-qual-reassign.log"
wait
echo $files > SAMPLES
./copy-from-webdav.sh
wait
localfiles=`cat LOCAL_SAMPLES`
grep '^Reads in allinfo before quality filter' $localfiles > report-ALLINFO.txt
echo "Wrote report-ALLINFO.txt"
grep '^Reads in allinfo after quality filter' $localfiles > report-ALLINFO-FILTER.txt
echo "Wrote report-ALLINFO-FILTER.txt"
grep '^Nr of clones' $localfiles > report-CLONES.txt
echo "Wrote report-CLONES.txt"
grep '^Nr of dominant clones' $localfiles > report-DOMINANT.txt
echo "Wrote report-DOMINANT.txt"
wait
rm -f $localfiles

# First argument: were additional mids used, yes or no
python2 report-combine-all.py ${BARCODES} ${INFO}
wait

echo "FINISHED"
