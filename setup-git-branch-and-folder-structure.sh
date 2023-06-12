#!/bin/bash

# Description: use this script to automagically create a branch for this specific dataset
#              and to create the FSS on the webdav server at SurfSara

# Fill in the ENV.sh, then run this script

source ENV.sh

# Checkout latest FSS structure
rm -rf Reproducibility
git clone git@github.com:EDS-Bioinformatics-Laboratory/Reproducibility.git

# Rename project directory name and analysis directory name
mv Reproducibility/_LATEST-ENCORE-TEMPLATE/yyyymmdd_ProjectName/ ${RUN}
mv ${RUN}/Processing/yyyymmdd_NameOfDataAnalysis/ ${RUN}/Processing/$OUTDIR/

# Create directory if it is not there yet, the "raw", "reports" and "final" directories are made via execute-all.sh
## LEGE DIRECTORIES WORDEN NIET NAAR WEBDAV GEKOPIEERD MET RCLONE
## TO DO: deze aanmaken in het execute-all.sh script
#mkdir -p $RUN/Processing/$OUTDIR/Settings # was: parameters
#mkdir -p $RUN/Processing/$OUTDIR/Results/clones
#mkdir -p $RUN/Processing/$OUTDIR/Results/run-report
#mkdir -p $RUN/Processing/$OUTDIR/Results/similarity
#mkdir -p $RUN/Processing/$OUTDIR/Results/shared-clones

# # Create a branch for this specific dataset
mybranch=${RUN}
#git branch ${mybranch}
#git checkout ${mybranch}

# Create a github.md file
echo "# Git" > $RUN/Processing/github.md
echo "" >> $RUN/Processing/github.md
echo "\`\`\`" >> $RUN/Processing/github.md
echo "git clone git@github.com:EDS-Bioinformatics-Laboratory/reseda.git" >> $RUN/Processing/github.md
echo "cd reseda" >> $RUN/Processing/github.md
echo "git branch ${mybranch}" >> $RUN/Processing/github.md
echo "git checkout ${mybranch}" >> $RUN/Processing/github.md
echo "git pull origin ${mybranch}" >> $RUN/Processing/github.md
echo "\`\`\`" >> $RUN/Processing/github.md

# Create a README.md file in the Processing/yyyymmdd_NameOfDataAnalysis directory
echo "## Code and settings" > $RUN/Processing/$OUTDIR/README.md
echo "" >> $RUN/Processing/$OUTDIR/README.md
echo "* [Settings/](Settings/) - files that are used as input for RESEDA" >> $RUN/Processing/$OUTDIR/README.md
echo "* [Code/](Code/) - information about code and versions" >> $RUN/Processing/$OUTDIR/README.md
echo "" >> $RUN/Processing/$OUTDIR/README.md
echo "## Standard output from RESEDA" >> $RUN/Processing/$OUTDIR/README.md
echo "" >> $RUN/Processing/$OUTDIR/README.md
echo "* [final/](final/) - contains all_info (per read) and clones files per sample" >> $RUN/Processing/$OUTDIR/README.md
echo "* [raw/](raw/) - intermediate files" >> $RUN/Processing/$OUTDIR/README.md
echo "* [reports/](reports/) - raw report files per sample" >> $RUN/Processing/$OUTDIR/README.md
echo "" >> $RUN/Processing/$OUTDIR/README.md
echo "## Post analysis with jupyter notebooks" >> $RUN/Processing/$OUTDIR/README.md
echo "" >> $RUN/Processing/$OUTDIR/README.md
echo "See git branch ${mybranch}" >> $RUN/Processing/$OUTDIR/README.md
echo "" >> $RUN/Processing/$OUTDIR/README.md
echo "* [run-report/](run-report/) - Performance of this run (e.g. identified MIDs, CDR3's, etc)" >> $RUN/Processing/$OUTDIR/README.md
echo "* [similarity/](similarity/) - Bray-Curtis analysis" >> $RUN/Processing/$OUTDIR/README.md
echo "* [shared-clones/](shared-clones/) - Shared clones analysis" >> $RUN/Processing/$OUTDIR/README.md

# Stub for lab journal LabJournal.md
echo "# Data analysis ${RUN} with RESEDA" > $RUN/ProjectDocumentation/LabJournal.md
echo "" >> $RUN/ProjectDocumentation/LabJournal.md
echo "Git repository: reseda" >> $RUN/ProjectDocumentation/LabJournal.md
echo "" >> $RUN/ProjectDocumentation/LabJournal.md
echo "Git branch: ${mybranch} (see _GIT.md)" >> $RUN/ProjectDocumentation/LabJournal.md
echo "" >> $RUN/ProjectDocumentation/LabJournal.md
echo "## Data analysis" >> $RUN/ProjectDocumentation/LabJournal.md
echo "" >> $RUN/ProjectDocumentation/LabJournal.md
echo "According to the documentation in the git repository: [./DOCS/_build/html/](https://github.com/EDS-Bioinformatics-Laboratory/reseda/tree/master/DOCS/_build/html)" >> $RUN/ProjectDocumentation/LabJournal.md
echo "" >> $RUN/ProjectDocumentation/LabJournal.md
echo "## Post analysis" >> $RUN/ProjectDocumentation/LabJournal.md
echo "" >> $RUN/ProjectDocumentation/LabJournal.md
echo "\`\`\`" >> $RUN/ProjectDocumentation/LabJournal.md
echo "./report-ALL.sh PARAMETERS" >> $RUN/ProjectDocumentation/LabJournal.md
echo "python ConcatenateCloneFilesBatch.py PARAMETERS" >> $RUN/ProjectDocumentation/LabJournal.md
echo "\`\`\`" >> $RUN/ProjectDocumentation/LabJournal.md
echo "" >> $RUN/ProjectDocumentation/LabJournal.md
echo "## Notebooks" >> $RUN/ProjectDocumentation/LabJournal.md
echo "" >> $RUN/ProjectDocumentation/LabJournal.md
echo "\`\`\`" >> $RUN/ProjectDocumentation/LabJournal.md
echo "NOTEBOOKS/SampleSimilarity.ipynb" >> $RUN/ProjectDocumentation/LabJournal.md
echo "NOTEBOOKS/SharedClonesDirection.ipynb" >> $RUN/ProjectDocumentation/LabJournal.md
echo "\`\`\`" >> $RUN/ProjectDocumentation/LabJournal.md

# # Upload files
cp _README-code.md $RUN/Processing/$OUTDIR/CodeDocumentation/README-code.md
rclone copy $RUN remote:amc-immunogenomics/RUNS/$RUN

# Some comment for the one who runs this script
echo "Created directory ${RUN} on the webdav server."
echo ""
echo "Please upload parameter files (PiCaS tokens) to ${WEBDAV}/Processing/${OUTDIR}/Settings/"
echo "Upload clones (vjcdr3-clones* cdr3-clones* assign-*, made with ConcatenateCloneFilesBatch.py) to ${WEBDAV}/Processing/${OUTDIR}/Results/clones/"
echo "Upload run reports (made by report-ALL.sh) to ${WEBDAV}/Processing/${OUTDIR}/Results/run-report/"
echo "Upload Bray-Curtis report (made with NOTEBOOKS/SampleSimilarity.ipynb) to ${WEBDAV}/Processing/${OUTDIR}/Results/similarity/"
echo "Upload shared-clones report (made with NOTEBOOKS/SharedClonesDirection.ipynb) to ${WEBDAV}/Processing/${OUTDIR}/Results/shared-clones/"
echo ""
echo "You are now in git branch ${mybranch}"
echo "Do not forget to commit and push all the dataset specific changes to github!!"

