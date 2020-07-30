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
    echo "    -c --cregion         default: no (sequences contain the C-region, align and mask this region)"
    echo "    -p --protocol        single|paired, default: paired"
    echo "    -o --outdir          default: results-tbcell"
    echo "    -b --barcodes        yes|no, were extra internal barcodes used? default=yes"
    echo "    -u --umis            yes|no, does sequence contain UMI? default=yes"
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
        -c|--cregion)
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
    runcmd ./copy-from-beehub.sh
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
    runcmd ./run-fastqc.sh ${samples}
    wait

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
runcmd ./run-fastqc.sh ${samples}
wait

# Align against the C-region and mask the sequences
if [ "${CREGION}" == "yes" ]; then
    runcmd ./batch-align.sh IGHC_CH12_human.fasta ${samples}
    samfiles=`ls split/*-IGHC_CH12_human.sam`
    python MaskSequences.py ${samfiles}
    samfiles=`ls split/*.masked.sam`
    for samfile in ${samfiles}; do
        java -Djava.io.tmpdir=./tmp -jar picard-tools-1.126/picard.jar SamToFastq I=${samfile} F=${samfile}.fastq
        gzip ${samfile}
    done
    samples=`ls split/*.masked.sam.fastq.gz`
fi

# Extract the CDR3 sequence
set_status ${ip} "RUNNING" "${CELLTYPE} Extracting CDR3's"
runcmd python2 TranslateAndExtractCdr3.py -c ${CELLTYPE} -m ${MISMATCHES} ${samples}
wait

# Count lines of CDR3.csv files
set_status ${ip} "RUNNING" "${CELLTYPE} Select correct MIDs"
runcmd wc -l split/*${CELLTYPE}-CDR3.csv > wc-${ip}.txt
wait
runcmd python2 select-correct-mids.py ${BARCODES} wc-${ip}.txt > mv-samples-with-correct-mid.sh
wait
mkdir split/correct-mid
wait
cd split
runcmd bash ../mv-samples-with-correct-mid.sh
wait
mv *.assembled-report.txt correct-mid
cd ..

samples=`ls split/correct-mid/*.fastq.gz`

# Align sequences against IMGT and call SNPs
set_status ${ip} "RUNNING" "${CELLTYPE} Aligning sequences"
for ref in $refs; do
   runcmd ./batch-align.sh ${ref} ${samples}
done
wait

# Get SNPs from SAM files
set_status ${ip} "RUNNING" "Determine SNPs from the SAM files" # creates file: ${prefix}-${refprefix}-e-clean.sam.mut.txt
runcmd python3 MutationsFromSam.py *-e-clean.sam
wait

### Generate reports ###

mkdir final

# For each sample; do
set_status ${ip} "RUNNING" "${CELLTYPE} Combining results"
for sample in ${samples}; do
    mydir=`dirname ${sample}`
    prefix=`basename ${sample} .fastq.gz`

    # Combine MID, CDR3, V, J and sequence information
    midFile=`echo ${mydir}/${prefix}|perl -ne 's/(.+)-.+$/$1-report.txt/;print;'`
    cdr3File=${sample}-${CELLTYPE}-CDR3.csv
    vFile=${prefix}-${v}-easy-import.txt
    jFile=${prefix}-${j}-easy-import.txt
    seqFile=${sample}-${CELLTYPE}.csv
    extraFile=${sample}-${CELLTYPE}-extra.txt
    allinfoFile="final/${prefix}-${CELLTYPE}-all_info.csv"
    cloneFile="final/${prefix}-${CELLTYPE}-clones.csv"
    cloneSubsFile="final/${prefix}-${CELLTYPE}-clones-subs.csv"
    cloneMainsFile="final/${prefix}-${CELLTYPE}-clones-mains.csv"
    totalFile="final/${prefix}-${CELLTYPE}-productive.txt"
    echo "### runcmd python2 combine-immuno-data.py ${midFile} ${cdr3File} ${vFile} ${jFile} ${seqFile} ${extraFile} ${allinfoFile} ${cloneFile} ${cloneSubsFile} ${cloneMainsFile} ${totalFile} ###"
    runcmd python2 combine-immuno-data.py ${midFile} ${cdr3File} ${vFile} ${jFile} ${seqFile} ${extraFile} ${allinfoFile} ${cloneFile} ${cloneSubsFile} ${cloneMainsFile} ${totalFile}
    wait

    # Integrate allinfo file with V and J mutation information, if it fails it will just continue with the next sample, creates a clones file
    vMutFile=${prefix}-${v}-e-clean.sam.mut.txt
    jMutFile=${prefix}-${j}-e-clean.sam.mut.txt
    python3 MutationAnalysisVJ.py -a ${allinfoFile} -v ${vMutFile} -j ${jMutFile}
    wait

    # Reassign V genes based on the created clones file above, creates a new clones file
    cloneMutFile=final/${prefix}-${CELLTYPE}-clones-mut-sites.csv # this is the result of the script MutationAnalysisVJ.py
    python3 ReassignGenes.py -c ${cloneMutFile} -a ${allinfoFile} # creates a file with extension -clones-mut-sites-reassigned.csv and -allinfo-filtered.csv
done

# Move results to 'final'
mv split/correct-mid/* final
mv final/*-report.txt split
wait

# Correct V gene assignments (OLD, can be removed when new procedure is correct)
set_status ${ip} "RUNNING" "${CELLTYPE} Re-assign V genes"
runcmd python2 ReAssignVGenes.py final/*-all_info.csv
wait
# Move files to final
mv *.rr.* final
wait

# Do mutation analysis if reference is IGH_HUMAN (OLD, can be removed when new procedure is correct)
if [[ ${CELLTYPE} -eq "IGH_HUMAN" ]]; then
    set_status ${ip} "RUNNING" "${CELLTYPE} Mutation analysis"
    samples=`ls final/*-IGH_HUMAN-all_info.csv.rr.all_info.csv`
    for sample in ${samples}; do
        prefix=`basename ${sample} -IGH_HUMAN-all_info.csv.rr.all_info.csv`
        runcmd Rscript MutationAnalysisVJ.R indir=\'.\' outdir=\'final\' V.file=\'${prefix}-IGHV_human-e-clean.sam.mut.txt\' J.file=\'${prefix}-IGHJ_human-e-clean.sam.mut.txt\' CDR3.file=\'final/${prefix}-IGH_HUMAN-all_info.csv.rr.all_info.csv\'
    done
fi

# Make output directories
mkdir ${beehub_mount}
mkdir ${beehub_mount}/${RESULTSDIR}
mkdir ${beehub_mount}/${RESULTSDIR}/raw
#mkdir ${beehub_mount}/${RESULTSDIR}/raw/correct-mid
mkdir ${beehub_mount}/${RESULTSDIR}/reports
mkdir ${beehub_mount}/${RESULTSDIR}/final
wait

# Transfer data to Beehub
set_status ${ip} "RUNNING" "Transferring ${CELLTYPE} data to Webdav server"
runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/reports/ *-pear.log *-pear.err *.quality-filter.log wc-*.txt versions-*
runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/reports/ split/*.primers.count.txt split/*-report.txt split/*-midcount.txt split/*-extra.txt
runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/reports/ final/*-productive.txt final/*.log

#runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/raw/ *.sam *.snp.csv *.mut.txt *.short*.assembled.fastq.gz
#runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/raw/ split/*.fastq.gz split/*_fastqc.zip split/*-alt-V-CDR3.csv split/*-alt-J-CDR3.csv
#runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/raw/correct-mid/ final/*L001* final/*mutations*
runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/raw/ final/*CDR3*.csv final/*discarded*.txt

runcmd ./copy-to-webdav.sh ${beehub_web}/${RESULTSDIR}/final/ final/*.rr.* final/*mutations* final/*-allinfo-filtered.csv final/*-clones-mut-sites.csv final/*-clones-mut-sites-reassigned.csv

wait

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"

# Set status when analysis is finished
set_status ${ip} "FINISHED" "${CELLTYPE} finished in ${difftime} seconds"

exit 0
