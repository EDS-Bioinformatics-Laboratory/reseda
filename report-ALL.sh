#!/bin/bash

################## Get arguments #################

function show_help {
    echo "USAGE: ./report-ALL.sh [options]"
    echo "    -r --run             required: runNN-YYYYMMDD-miseq"
    echo "    -i --info            required: yyyymmdd_RUNnn_Datasheet-new.json"
    echo "    -b --barcodes        default: yes (other option: no)"
    echo "    -o --outdir          default: results-tbcell"
    echo "    -p --prefix          name of the rclone remote, default: remote"
    exit
}

OUTDIR="results-tbcell"
BARCODES="yes"
PREFIX="remote"

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
        -p|--prefix)
        PREFIX="$2"
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
echo PREFIX          = "${PREFIX}"
# echo POSITIONAL     = "${POSITIONAL}"
# echo rest           = "$@"

################## End get arguments #############

remotepath="${PREFIX}:RUNS/${RUN}/Processing/${OUTDIR}/Results"
tmpdir="report-tmp"
echo "Remote path: ${remotepath}"

# Download files from the webdav server with rclone
# Usage: fetch_files subdir pattern [pattern ...]   e.g.: fetch_files reports '*-pear.log'
# Sets: localfiles (paths of the downloaded files), returns 1 if no files were found
function fetch_files {
    local subdir=$1
    shift
    local includes=()
    for pattern in "$@"; do
        includes+=(--include "${pattern}")
    done

    rm -rf ${tmpdir}
    mkdir -p ${tmpdir}
    echo "rclone copy ${remotepath}/${subdir} ${tmpdir} ${includes[*]}"
    rclone copy ${remotepath}/${subdir} ${tmpdir} "${includes[@]}" --progress --transfers=4 --checkers=8 --retries=3
    wait
    localfiles=`ls ${tmpdir}/* 2> /dev/null`
    if [ -z "${localfiles}" ]; then
        echo "WARNING: no files found in ${remotepath}/${subdir} for: $*"
        return 1
    fi
    return 0
}

# grep -H: always print the file name, also when there is only one file (the report-combine-all.py parser needs it)

if fetch_files reports '*-pear.log'; then
    grep -H '^Assembled reads \.' $localfiles > report-PEAR.txt
else
    > report-PEAR.txt
fi
wait
echo "Wrote report-PEAR.txt"

# Same format as the output of "more file1 file2 ...": a header with the file name per file
if fetch_files reports '*-midcount.txt'; then
    for file in $localfiles; do
        echo "::::::::::::::"
        echo ${file}
        echo "::::::::::::::"
        cat ${file}
    done > report-MIDs.txt
else
    > report-MIDs.txt
fi
wait
echo "Wrote report-MIDs.txt"

if fetch_files reports '*HUMAN-report.txt' '*MOUSE-report.txt' '*RHESUS-report.txt'; then
    grep -H '^4' $localfiles > report-CDR3.txt
    grep -H '^5' $localfiles > report-ALT-V.txt
    grep -H '^6' $localfiles > report-ALT-J.txt
else
    > report-CDR3.txt
    > report-ALT-V.txt
    > report-ALT-J.txt
fi
wait
echo "Wrote report-CDR3.txt"
echo "Wrote report-ALT-V.txt"
echo "Wrote report-ALT-J.txt"

if fetch_files reports '*-productive.txt'; then
    grep -H 'Unique reads with V and J in all_info' $localfiles > report-PRODUCTIVE.txt
else
    > report-PRODUCTIVE.txt
fi
wait
echo "Wrote report-PRODUCTIVE.txt"

if fetch_files raw '*.sam'; then
    python3 report-alignments.py ${localfiles}
    wait
    echo "Wrote report-ALIGNED-*"
fi

if fetch_files final '*.rr.clones_subs.csv'; then
    python3 report-after-v-reassignment.py $localfiles
else
    > report-AFTER-V-REASSIGNMENT.txt
fi
wait
echo "Wrote report-AFTER-V-REASSIGNMENT.txt"

if fetch_files reports '*.quality-filter.log'; then
    perl -ne 'print if ! m/^datafile/;' $localfiles > report-QUALITY-FILTER.txt
else
    > report-QUALITY-FILTER.txt
fi
wait
echo "Wrote report-QUALITY-FILTER.txt"

if fetch_files reports '*-qual-reassign.log'; then
    grep -H '^Reads in allinfo before quality filter' $localfiles > report-ALLINFO.txt
    grep -H '^Reads in allinfo after quality filter' $localfiles > report-ALLINFO-FILTER.txt
    grep -H '^Nr of clones' $localfiles > report-CLONES.txt
    grep -H '^Nr of dominant clones' $localfiles > report-DOMINANT.txt
else
    > report-ALLINFO.txt
    > report-ALLINFO-FILTER.txt
    > report-CLONES.txt
    > report-DOMINANT.txt
fi
wait
echo "Wrote report-ALLINFO.txt"
echo "Wrote report-ALLINFO-FILTER.txt"
echo "Wrote report-CLONES.txt"
echo "Wrote report-DOMINANT.txt"

rm -rf ${tmpdir}

# First argument: were additional mids used, yes or no
python3 report-combine-all.py ${BARCODES} ${INFO}
wait

echo "FINISHED"
