#!/bin/bash

# Load environment variables: RUN OUTDIR WEBDAV
source ENV.sh
wait

# Create folder structure

# Mount the Beehub webdav server and configure the location
webdav_mount="/mnt/immunogenomics/RUNS/${RUN}"
webdav_url="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/${RUN}"

mkdir -p ${webdav_mount}/Data/Meta
mkdir -p ${webdav_mount}/Data/Raw
mkdir -p ${webdav_mount}/Data/Preprocessed

mkdir -p ${webdav_mount}/Processing/${OUTDIR}/Documentation     # upload versions* files here
mkdir -p ${webdav_mount}/Processing/${OUTDIR}/ResultsPerSample  # what is now uploaded to "final"
mkdir -p ${webdav_mount}/Processing/${OUTDIR}/ResultsRun        # upload cdr3-clones and report-all* here
mkdir -p ${webdav_mount}/Processing/${OUTDIR}/LogsAndIntermediateFiles  # what is now uploaded to "raw"
mkdir -p ${webdav_mount}/Processing/${OUTDIR}/Settings          # tokens
mkdir -p ${webdav_mount}/Processing/${OUTDIR}/Notebooks         # subdirs, such as "similarity", "shared-clones", "rarefaction"

# Create/Upload README.md and CONTENT.md files
touch README-analysis.md    # Lab journal
touch GIT-master.md         # Link to code and instructions on how to download
touch GIT-data-analysis.md  # Git branch for this specific analysis

## What to change in manual steps (update the documentation/README.md!):

# Upload raw data and meta data:
# - copy-basespace-data-to-beehub.py (to Data/Raw)
# - manual: upload Datasheet to Data/Meta
# - when datasheet is altered, files are renamed, files are converted: put it in Data/Preprocessed including the code or reference to the code

## What to change in the code?

# execute-all.sh:
# - remove creation of researchdrive directories
# - upload results according to the data structure above

# MakeSamplesFiles.py:
# - change the instruction in the help text
# - ToposCreateTokens.py??
# - picas/createTokens.py??
# - pilot.py??

