#!/bin/bash

# Description: use this script to automagically create a branch for this specific dataset
#              and upload general README, CONTENT, and GIT files to the webdav server

# Mount the research drive under /mnt/immunogenomics/
# Fill in the ENV.sh, then run this script

source ENV.sh

# Create directory if it is not there yet, the "raw", "reports" and "final" directories are made via execute-all.sh
mkdir -p /mnt/immunogenomics/RUNS/$RUN/data/raw
mkdir -p /mnt/immunogenomics/RUNS/$RUN/data/meta
mkdir -p /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/code
mkdir -p /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/parameters
mkdir -p /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/run-report
mkdir -p /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/similarity
mkdir -p /mnt/immunogenomics/RUNS/$RUN/$OUTDIR/shared-clones

# # Create a branch for this specific dataset
mybranch=${RUN}
git branch ${mybranch}
git checkout ${mybranch}

# Create a _GIT.md file
echo "# Git" > _GIT.md
echo "" >> _GIT.md
echo "\`\`\`" >> _GIT.md
echo "git clone git@github.com:EDS-Bioinformatics-Laboratory/reseda.git" >> _GIT.md
echo "cd reseda" >> _GIT.md
echo "git branch ${mybranch}" >> _GIT.md
echo "git checkout ${mybranch}" >> _GIT.md
echo "git pull origin ${mybranch}" >> _GIT.md
echo "\`\`\`" >> _GIT.md

# Create a _CONTENT.md file
echo "## Code and settings" > _CONTENT.md
echo "" >> _CONTENT.md
echo "* [parameters/](parameters/) - files that are used as input for RESEDA" >> _CONTENT.md
echo "* [code/](code/) - information about code and versions" >> _CONTENT.md
echo "" >> _CONTENT.md
echo "## Standard output from RESEDA" >> _CONTENT.md
echo "" >> _CONTENT.md
echo "* [final/](final/) - contains all_info (per read) and clones files per sample" >> _CONTENT.md
echo "* [raw/](raw/) - intermediate files" >> _CONTENT.md
echo "* [reports/](reports/) - raw report files per sample" >> _CONTENT.md
echo "" >> _CONTENT.md
echo "## Post analysis with jupyter notebooks" >> _CONTENT.md
echo "" >> _CONTENT.md
echo "See git branch ${mybranch}" >> _CONTENT.md
echo "" >> _CONTENT.md
echo "* [run-report/](run-report/) - Performance of this run (e.g. identified MIDs, CDR3's, etc)" >> _CONTENT.md
echo "* [similarity/](similarity/) - Bray-Curtis analysis" >> _CONTENT.md
echo "* [shared-clones/](shared-clones/) - Shared clones analysis" >> _CONTENT.md

# Stub for lab journal _LabJournal.md
echo "# Data analysis ${RUN} with RESEDA" > _LabJournal.md
echo "" >> _LabJournal.md
echo "Git repository: reseda" >> _LabJournal.md
echo "" >> _LabJournal.md
echo "Git branch: ${mybranch} (see _GIT.md)" >> _LabJournal.md
echo "" >> _LabJournal.md
echo "## Data analysis" >> _LabJournal.md
echo "" >> _LabJournal.md
echo "According to the documentation in the git repository: [./DOCS/_build/html/](https://github.com/EDS-Bioinformatics-Laboratory/reseda/tree/master/DOCS/_build/html)" >> _LabJournal.md
echo "" >> _LabJournal.md
echo "## Post analysis" >> _LabJournal.md
echo "" >> _LabJournal.md
echo "\`\`\`" >> _LabJournal.md
echo "./report-ALL.sh PARAMETERS" >> _LabJournal.md
echo "python ConcatenateCloneFilesBatch.py PARAMETERS" >> _LabJournal.md
echo "\`\`\`" >> _LabJournal.md
echo "" >> _LabJournal.md
echo "## Notebooks" >> _LabJournal.md
echo "" >> _LabJournal.md
echo "\`\`\`" >> _LabJournal.md
echo "NOTEBOOKS/SampleSimilarity.ipynb" >> _LabJournal.md
echo "NOTEBOOKS/SharedClonesDirection.ipynb" >> _LabJournal.md
echo "\`\`\`" >> _LabJournal.md

# # Upload files
./copy-to-webdav.sh ${WEBDAV}/${OUTDIR}/ _GIT.md _CONTENT.md _LabJournal.md
./copy-to-webdav.sh ${WEBDAV}/${OUTDIR}/code/ _README-code.md

# Some comment for the one who runs this script
echo "Created directory ${RUN}/${OUTDIR} on the webdav server."
echo "Uploaded _GIT.md _CONTENT.md _LabJournal.md to the webdav server"
echo ""
echo "Please upload parameter files (PiCaS tokens) to ${WEBDAV}/${OUTDIR}/parameters/"
echo "Upload clones (vjcdr3-clones* cdr3-clones* assign-*, made with ConcatenateCloneFilesBatch.py) to ${WEBDAV}/${OUTDIR}/"
echo "Upload run reports (made by report-ALL.sh) to ${WEBDAV}/${OUTDIR}/run-report/"
echo "Upload Bray-Curtis report (made with NOTEBOOKS/SampleSimilarity.ipynb) to ${WEBDAV}/${OUTDIR}/similarity/"
echo "Upload shared-clones report (made with NOTEBOOKS/SharedClonesDirection.ipynb) to ${WEBDAV}/${OUTDIR}/shared-clones/"
echo ""
echo "You are now in git branch ${mybranch}"
echo "Do not forget to commit and push all the dataset specific changes to github!!"
