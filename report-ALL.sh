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

files="/mnt/immunogenomics/RUNS/${RUN}/${OUTDIR}/reports/*.log"
grep '^Assembled reads \.' $files > report-PEAR.txt
echo "Wrote report-PEAR.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/${OUTDIR}/reports/*-midcount.txt"
more $files > report-MIDs.txt
echo "Wrote report-MIDs.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/${OUTDIR}/reports/*HUMAN-report.txt /mnt/immunogenomics/RUNS/${RUN}/${OUTDIR}/reports/*MOUSE-report.txt"
grep '^4' $files > report-CDR3.txt
echo "Wrote report-CDR3.txt"
grep '^5' $files > report-ALT-V.txt
echo "Wrote report-ALT-V.txt"
grep '^6' $files > report-ALT-J.txt
echo "Wrote report-ALT-J.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/${OUTDIR}/reports/*-productive.txt"
grep 'Unique reads with V and J in all_info' $files > report-PRODUCTIVE.txt
echo "Wrote report-PRODUCTIVE.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/${OUTDIR}/raw/*.sam"
python report-alignments.py ${files}
echo "Wrote report-ALIGNED-*"

python2 report-after-v-reassignment.py /mnt/immunogenomics/RUNS/${RUN}/${OUTDIR}/final/*.rr.clones_subs.csv
echo "Wrote report-AFTER-V-REASSIGNMENT.txt"

files="/mnt/immunogenomics/RUNS/${RUN}/${OUTDIR}/reports/*.quality-filter.log"
perl -ne 'print if ! m/^datafile/;' $files > report-QUALITY-FILTER.txt
echo "Wrote report-QUALITY-FILTER.txt"

# First argument: were additional mids used, yes or no
python2 report-combine-all.py ${BARCODES} ${INFO}
wait

echo "FINISHED"
