#!/bin/bash

# Preparation:
# mv reference/* .
# mv hlaforest/scripts/* .

# Create a file with all the fastq files:
# ls TESTDATA/* > SAMPLES

################## Get arguments #################

function show_help {
    echo "USAGE: ./execute-all.sh [options]"
    echo "    -r --run             required: runNN-YYYYMMDD-miseq"
    echo "    -l --location        local|webdav, default: webdav"
    echo "    -m --mids            default: MIDS-miseq-umi.txt (older version is MIDS-miseq.txt)"
    echo "    -org --organism      human|mouse, default: human"
    echo "    -cell --cell         IGH|IGK|IGL|TRA|TRB, default:IGH"
    echo "    -celltype --celltype IGH_HUMAN|TRB_MOUSE|etc, default: IGH_HUMAN"
    echo "    -mm --mismatches     default: 0 (mismatches allowed in CDR3 motifs)"
    echo "    -cregion --cregion         default: no (sequences contain the C-region, align and mask this region)"
    echo "    -p --protocol        single|paired, default: paired"
    echo "    -o --outdir          default: results-tbcell"
    echo "    -b --barcodes        yes|no, were extra internal barcodes used? default=yes"
    echo "    -u --umis            yes|roche|race|no, does sequence contain UMI? default=yes"
    exit
}

LOCATION="webdav"
MIDS="MIDS-miseq-umi.txt"
ORGANISM="human"
CELL="IGH"
CELLTYPE="IGH_HUMAN"
MISMATCHES=0
CREGION="no"
PROTOCOL="paired"
RESULTSDIR="results-tbcell"
BARCODES="yes"
UMIS="yes"

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
        -l|--location)
        LOCATION="$2"
        shift # past argument
        shift # past value
        ;;
        -m|--mids)
        MIDS="$2"
        shift # past argument
        shift # past value
        ;;
        -org|--organism)
        ORGANISM="$2"
        shift # past argument
        shift # past value
        ;;
        -cell|--cell)
        CELL="$2"
        shift # past argument
        shift # past value
        ;;
        -celltype|--celltype)
        CELLTYPE="$2"
        shift # past argument
        shift # past value
        ;;
        -mm|--mismatches)
        MISMATCHES="$2"
        shift # past argument
        shift # past value
        ;;
        -cregion|--cregion)
        CREGION="$2"
        shift # past argument
        shift # past value
        ;;
        -p|--protocol)
        PROTOCOL="$2"
        shift # past argument
        shift # past value
        ;;
        -o|--outdir)
        RESULTSDIR="$2"
        shift # past argument
        shift # past value
        ;;
        -b|--barcodes)
        BARCODES="$2"
        shift # past argument
        shift # past value
        ;;
        -u|--umis)
        UMIS="$2"
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

echo RUN             = "${RUN}"
echo LOCATION        = "${LOCATION}"
echo MIDS            = "${MIDS}"
echo ORGANISM        = "${ORGANISM}"
echo CELL            = "${CELL}"
echo CELLTYPE        = "${CELLTYPE}"
echo MISMATCHES      = "${MISMATCHES}"
echo CREGION         = "${CREGION}"
echo PROTOCOL        = "${PROTOCOL}"
echo RESULTSDIR      = "${RESULTSDIR}"
echo BARCODES        = "${BARCODES}"
echo UMIS            = "${UMIS}"
# echo POSITIONAL      = "${POSITIONAL}"
# echo rest            = "$@"

################## End get arguments #############

# Reference sequences
refs="${CELL}V_${ORGANISM}.fasta ${CELL}J_${ORGANISM}.fasta"
v="${CELL}V_${ORGANISM}"
j="${CELL}J_${ORGANISM}"

# Mount the Beehub webdav server and configure the location
beehub_mount="/mnt/immunogenomics/RUNS/${RUN}"
beehub_web="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/${RUN}"

# Then run ./execute-all.sh


#######################

starttime=`date +%s`
ip_address=`hostname -I`
ips=($ip_address)
ip=${ips[0]}

thisdir=`pwd`

function runcmd {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "ERROR ${CELLTYPE} with $1" >&2
        set_status ${ip} "ERROR" "Error ${CELLTYPE} with $1"
        exit 1
    fi
    return $status
}

# Job monitoring. To use this you need https://bitbucket.org/barbera/progress
function set_status {
    local ip=$1
    local stat=$2
    local message=$3
    cd ../progress
    python3 set-status.py ip:"${ip}" status:"${stat}" message:"${message}"
    cd $thisdir
}

runcmd ./log-versions.sh

set_status ${ip} "RUNNING" "Started ${MIDS} ${CELLTYPE} analysis on ${starttime}"

### Change this if you work with local/remote files ###

if [ "${LOCATION}" == "webdav" ]; then
    # Remote files:
    runcmd ./copy-from-webdav.sh
    wait
else
    # Local files:
    cp SAMPLES LOCAL_SAMPLES
    wait
fi


#######################################################

samples=`cat LOCAL_SAMPLES`  # get all arguments
r1_samples=`grep R1_001 LOCAL_SAMPLES`

### Analysis on raw fastq files ###
if [ "${PROTOCOL}" != "single" ]; then
    # FastQC
    #runcmd ./run-fastqc.sh ${samples}
    #wait

    # Pairwise assembly
    set_status ${ip} "RUNNING" "${CELLTYPE} Pairwise assembly"
    runcmd ./batch-pear.sh ${r1_samples}
    wait

    samples=`ls *.assembled.fastq.gz`
fi

### Continue with assembled fastq files ###

# # Split on sequence length
# set_status ${ip} "RUNNING" "${CELLTYPE} Split sequences on length"
# runcmd python2 FastqSplitOnSequenceLength.py -l 270 ${samples}
# wait
#
# # New sample list
# samples=`cat SAMPLES_long`

# Split on MID
set_status ${ip} "RUNNING" "${CELLTYPE} Sorting sequences per MID"
runcmd python2 FastqSplitOnMid.py ${UMIS} ${MIDS} split ${samples}
wait

### Continue with the assembled, split per mid, fastq files ###

samples=`ls split/*.fastq.gz`

# FastQC report
#runcmd ./run-fastqc.sh ${samples}
#wait

# Align against the C-region and mask the sequences
if [ "${CREGION}" == "yes" ]; then
    runcmd ./batch-align.sh IGHC_CH12_human.fasta ${samples}
    wait
    mv *-IGHC_CH12_human.sam split/
    wait
    samfiles=`ls split/*-IGHC_CH12_human.sam`
    python MaskSequences.py ${samfiles}
    wait
    # copy the original fastq files to directory "orig" and convert to tab
    mkdir -p orig/correct-mid
    cp ${samples} orig/
    runcmd python2 SeqToFastq.py fastq2tab orig/*.fastq.gz
    wait
    mv *.tab.csv orig
    # Rename the masked files back to the original fastq file names
    ls split/*-IGHC_CH12_human.masked.fastq.gz |perl -ne 's/\n//; $masked = $_; s/-IGHC_CH12_human.masked//; $orig = $_; print "rename $masked $orig\n"; rename $masked, $orig;'
    wait
    samples=`ls split/*.fastq.gz`
fi


### Generate reports ###

# Transfer data to Beehub
set_status ${ip} "RUNNING" "Transferring ${CELLTYPE} data to Webdav server"

runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/raw/ *.sam *.snp.csv *.mut.txt *.short*.assembled.fastq.gz
runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/raw/ split/*.fastq.gz
runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/code/ versions-*

wait

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"

# Set status when analysis is finished
set_status ${ip} "FINISHED" "${CELLTYPE} finished in ${difftime} seconds"

exit 0
