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
    echo "    -s --seqlength      default: 0 (threshold on sequence length)"
    echo "    -cregion --cregion         default: no (sequences contain the C-region, align and mask this region)"
    echo "    -p --protocol        single|paired, default: paired"
    echo "    -o --outdir          default: results-tbcell"
    echo "    -b --barcodes        yes|no, were extra internal barcodes used? default=yes"
    echo "    -u --umis            yes|roche|race|vasaseq|no, does sequence contain UMI? default=yes"
    exit
}

LOCATION="webdav"
MIDS="MIDS-miseq-umi.txt"
ORGANISM="human"
CELL="IGH"
CELLTYPE="IGH_HUMAN"
MISMATCHES=0
SEQLENGTH=0
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
        -s|--seqlength)
        SEQLENGTH="$2"
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
echo SEQLENGTH       = "${SEQLENGTH}"
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
    python set-status.py ip:"${ip}" status:"${stat}" message:"${message}"
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

# Split on sequence length
if [ "${SEQLENGTH}" -gt "0" ]; then
    set_status ${ip} "RUNNING" "${CELLTYPE} Split sequences on length"
    runcmd python FastqSplitOnSequenceLength.py -l ${SEQLENGTH} ${samples}
    wait

    # New sample list
    samples=`cat SAMPLES_long`
fi

# Split on MID
set_status ${ip} "RUNNING" "${CELLTYPE} Sorting sequences per MID"
runcmd python FastqSplitOnMid.py ${UMIS} ${MIDS} split ${samples}
wait

### Continue with the assembled, split per mid, fastq files ###

samples=`ls split/*.fastq.gz`

# FastQC report
runcmd ./run-fastqc.sh ${samples}
wait

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
    runcmd python SeqToFastq.py fastq2tab orig/*.fastq.gz
    wait
    mv *.tab.csv orig
    # Rename the masked files back to the original fastq file names
    ls split/*-IGHC_CH12_human.masked.fastq.gz |perl -ne 's/\n//; $masked = $_; s/-IGHC_CH12_human.masked//; $orig = $_; print "rename $masked $orig\n"; rename $masked, $orig;'
    wait
    samples=`ls split/*.fastq.gz`
fi

# Extract the CDR3 sequence
set_status ${ip} "RUNNING" "${CELLTYPE} Extracting CDR3's"
runcmd python TranslateAndExtractCdr3.py -c ${CELLTYPE} -m ${MISMATCHES} ${samples}
wait

# Count lines of CDR3.csv files
set_status ${ip} "RUNNING" "${CELLTYPE} Select correct MIDs"
runcmd wc -l split/*${CELLTYPE}-CDR3.csv > wc-${ip}.txt
wait
runcmd python select-correct-mids.py ${BARCODES} wc-${ip}.txt > mv-samples-with-correct-mid.sh
wait
mkdir -p split/correct-mid
wait
cd split
runcmd bash ../mv-samples-with-correct-mid.sh
wait
mv *.assembled-report.txt correct-mid
cd ..
if [ "${CREGION}" == "yes" ]; then
    cd orig
    runcmd bash ../mv-samples-with-correct-mid.sh
    wait
    cd ..
fi

samples=`ls split/correct-mid/*.fastq.gz`

# Align sequences against IMGT and call SNPs
set_status ${ip} "RUNNING" "${CELLTYPE} Aligning sequences"
for ref in $refs; do
   runcmd ./batch-align.sh ${ref} ${samples}
done
wait

# Get SNPs from SAM files
set_status ${ip} "RUNNING" "Determine SNPs from the SAM files" # creates file: ${prefix}-${refprefix}-e-clean.sam.mut.txt
runcmd python MutationsFromSam.py *-e-clean.sam
wait

### Generate reports ###

mkdir -p final

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
    echo "### runcmd python combine-immuno-data.py -m ${midFile} -c ${cdr3File} -v ${vFile} -j ${jFile} -s ${seqFile} -e ${extraFile} -o ${allinfoFile} -ocf ${cloneFile} -ocs ${cloneSubsFile} -ocm ${cloneMainsFile} -t ${totalFile} ###"
    runcmd python combine-immuno-data.py -m ${midFile} -c ${cdr3File} -v ${vFile} -j ${jFile} -s ${seqFile} -e ${extraFile} -o ${allinfoFile} -ocf ${cloneFile} -ocs ${cloneSubsFile} -ocm ${cloneMainsFile} -t ${totalFile}
    wait

    # Integrate allinfo file with V and J mutation information, if it fails it will just continue with the next sample, creates a clones file
    vMutFile=${prefix}-${v}-e-clean.sam.mut.txt
    jMutFile=${prefix}-${j}-e-clean.sam.mut.txt
    python MutationAnalysisVJ.py -a ${allinfoFile} -v ${vMutFile} -j ${jMutFile}
    wait

    # Reassign V genes based on the created clones file above, creates a new clones file
    cloneMutFile=final/${prefix}-${CELLTYPE}-clones-mut-sites.csv # this is the result of the script MutationAnalysisVJ.py
    python ReassignGenes.py -c ${cloneMutFile} -a ${allinfoFile} # creates a file with extension -clones-mut-sites-reassigned.csv and -allinfo-filtered.csv
done

# Move results to 'final'
mv split/correct-mid/* final
mv final/*-report.txt split
wait

# Correct V gene assignments (OLD, can be removed when new procedure is correct)
set_status ${ip} "RUNNING" "${CELLTYPE} Re-assign V genes"
runcmd python ReAssignVGenes.py final/*-all_info.csv
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
mkdir -p ${beehub_mount}/Processing/${RESULTSDIR}/Results/versions
mkdir -p ${beehub_mount}/Processing/${RESULTSDIR}/Results/raw
mkdir -p ${beehub_mount}/Processing/${RESULTSDIR}/Results/reports
mkdir -p ${beehub_mount}/Processing/${RESULTSDIR}/Results/final

# Also make these directories. The data analysis results for the directory below are uploaded manually
mkdir -p ${beehub_mount}/Processing/${RESULTSDIR}/Results/clones
mkdir -p ${beehub_mount}/Processing/${RESULTSDIR}/Results/run-report
mkdir -p ${beehub_mount}/Processing/${RESULTSDIR}/Results/similarity
mkdir -p ${beehub_mount}/Processing/${RESULTSDIR}/Results/shared-clones

wait

# Transfer data to Beehub
set_status ${ip} "RUNNING" "Transferring ${CELLTYPE} data to Webdav server"
runcmd ./copy-to-webdav.sh ${beehub_web}/Processing/${RESULTSDIR}/Results/reports/ *-pear.log *-pear.err *.quality-filter.log wc-*.txt report-*.txt
runcmd ./copy-to-webdav.sh ${beehub_web}/Processing/${RESULTSDIR}/Results/reports/ split/*.primers.count.txt split/*-report.txt split/*-midcount.txt split/*-extra.txt
runcmd ./copy-to-webdav.sh ${beehub_web}/Processing/${RESULTSDIR}/Results/reports/ final/*-productive.txt final/*.log

#runcmd ./copy-to-webdav.sh ${beehub_web}/Processing/${RESULTSDIR}/Results/raw/ *.sam *.snp.csv *.mut.txt *.short*.assembled.fastq.gz
#runcmd ./copy-to-webdav.sh ${beehub_web}/Processing/${RESULTSDIR}/Results/raw/ split/*.fastq.gz split/*_fastqc.zip split/*-alt-V-CDR3.csv split/*-alt-J-CDR3.csv
#runcmd ./copy-to-webdav.sh ${beehub_web}/Processing/${RESULTSDIR}/Results/raw/correct-mid/ final/*L001* final/*mutations*
runcmd ./copy-to-webdav.sh ${beehub_web}/Processing/${RESULTSDIR}/Results/versions/ versions-*
runcmd ./copy-to-webdav.sh ${beehub_web}/Processing/${RESULTSDIR}/Results/raw/ final/*CDR3*.csv final/*discarded*.txt

runcmd ./copy-to-webdav.sh ${beehub_web}/Processing/${RESULTSDIR}/Results/final/ final/*.rr.* final/*mutations* final/*-allinfo-filtered.csv final/*-allinfo-filtered-mut.csv final/*-clones-mut-sites.csv final/*-clones-mut-sites-reassigned.csv

# Transfer the split fastq files that were converted to tab
if [ "${CREGION}" == "yes" ]; then
    mkdir -p ${beehub_mount}/Processing/${RESULTSDIR}/Results/fastq2tab
    runcmd ./copy-to-webdav.sh ${beehub_web}/Processing/${RESULTSDIR}/Results/fastq2tab/ orig/correct-mid/*.tab.csv
fi

wait

endtime=`date +%s`
difftime=`expr ${endtime} - ${starttime}`
echo "FINISHED WITH EXECUTE-ALL IN $difftime seconds"

# Set status when analysis is finished
set_status ${ip} "FINISHED" "${CELLTYPE} finished in ${difftime} seconds"

exit 0
