# RESEDA - REpertoire SEquencing Data Analysis

Data analysis workflow for T- and B-cell receptor repertoire sequencing.
The workflow identifies clones and their frequency from next generation sequencing of repertoires and includes steps for quality control and bias correction.

## Authors

Barbera DC van Schaik - b.d.vanschaik@amsterdamumc.nl

Prof. dr. Antoine HC van Kampen - a.h.vankampen@amsterdamumc.nl

## Workflow

![workflow](workflow.png)

## Required software

* [Samtools](http://www.htslib.org/)

The software packages below are included in this repository for convenience. Please visit the websites for more recent versions and information about the licenses.

* [PEAR](https://cme.h-its.org/exelixis/web/software/pear/doc.html)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [BWA](http://bio-bwa.sourceforge.net/)
* [VarScan](http://varscan.sourceforge.net/)
* [Picard](https://broadinstitute.github.io/picard/)

## Other requirements

The older scripts only work in Python 2 (see execute-all.sh)
* Python 2.7 and 3.5 (or higher)
    * future (print_function)
    * sys
    * os
    * subprocess
    * gzip
    * re
    * math
    * random
    * biopython
    * regex
    * sqlite3
    * matplotlib
    * numpy
    * scipy
    * json
    * pandas
    * shutil
    * argparse
    * csv
    * couchdb (when using the PiCaS pilotjob server)
* R
    * plyr
* Bash

## Settings webdav server

These scripts assume that the data is on the SurfSara ResearchDrive (webdav server) and that the drive is mounted in /mnt/immunogenomics
* copy-from-beehub.sh
* copy-basespace-data-to-beehub.py
* execute-all.sh
* report-ALL.sh
* ConcatenateCloneFilesBatch.py
* MakeSamplesFiles.py
* VerifyBasespaceCopy.py

sudo mount -t davfs -o uid=bioinfo,gid=bioinfo,rw https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics /mnt/immunogenomics

Create a .netrc file in your home directory:
* machine researchdrive.surfsara.nl
* login YOUR_RESEARCHDRIVE_USER
* password YOUR_PASSWORD

## How to run - Using only the code in this repository

Copy/move the relevant (or all) database files from directories "reference", "reftables" and "mids" to the root directory of the repository (same directory as the execute-all.sh script)

The input files (in fastq format) can be specified by putting the paths in the file SAMPLES. Running ./execute-all.sh without arguments shows you which parameters can be set.

Example: ./execute-all.sh -r mytestrun -l local -m MIDS-miseq.txt -org human -cell IGH -celltype IGH_HUMAN -u no

The results will be on the machine where you run this script, and if the webdav settings are configured well a new directory with all the files will appear on the ResearchDrive:
https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/mytestrun

## How to run - Using PiCaS for sending jobs

Note: this is what Barbera does for each sequence run

PiCaS: https://picas.surfsara.nl:6984/_utils/database.html?tbcellrep-bschaik-hpc

Ansible: [../ansible-playbooks/Reseda/README.md](../ansible-playbooks/Reseda/README.md)

### Setting up the FSS data structure and create a new git branch

Note: this step can run on a Linux laptop or on a virtual machine. In the latter case you need to install the software on the VM first (with Ansible)

Edit the ``ENV.sh`` file to specify the new directory on the ResearchDrive. Then run:

```
./setup-git-branch-and-folder-structure.sh
```

### Preparation ###

Note: these steps can run on a Linux laptop or on a virtual machine. In the latter case you need to install the software on the VM first (with Ansible)

* Mount basespace. Instructions are in basespace.txt
* Specify the run and the basespace sub-directories as argument to copy-basespace-data-to-beehub.py and run it. The file basespace-copy-data.sh and basespace-calc-checksum.sh will be created. Run these scripts to copy the data from basespace to the ResearchDrive and to calculate the SHA1 sums
* Convert the MiSeq sample sheet with MetaData.py (creates a json file)
    * An example of a datasheet can be found in the run directories on the ResearchDrive (data/meta/)
* Mount the ResearchDrive webdav server
* Add extra information to the json file with MakeSamplesFiles.py (this will also make the SAMPLE-* files)
* Sort and split the SAMPLE-* files with: ./SortAndSplit.sh SAMPLE-* It does the following:
    * Sorts the SAMPLE-* files: sort SAMPLE-blah > SAMPLE-blah.sort
    * Makes manageable jobs by splitting the SAMPLE-\*.sort files, e.g.: split -l 20 SAMPLES-run13-human-BCRh.sort SAMPLES-run13-human-BCRh-
* Delete existing jobs from the PiCaS server:
    * ``python deleteTokens.py Monitor/locked``
    * ``python deleteTokens.py Monitor/todo``
    * ``python deleteTokens.py Monitor/done``
* Create PiCaS jobs with ToposCreateTokens.py (run with the -h option to see the arguments)
* Upload PiCaS jobs with picas/createTokens.py JSON-FILES (json files that were created in the previous step)
* Create a "view" in the database with picas/createViews.py

### Starting the jobs ###

* Start virtual machines for the analysis (in the SurfSara cloud webinterface)
* Ansible scripts are used to install the machines and to start the data analysis
    * See: [../ansible-playbooks/Reseda/README.md](../ansible-playbooks/Reseda/README.md)

### When all jobs are finished ###

* In the jobs the data is automatically transferred to the ResearchDrive webdav server
* Execute ConcatenateCloneFilesBatch.py to generate a bash script for concatenating clone files per project+organism+cell_type (run the generated script)
* Run report-ALL.sh to generate reports about the sequence run (help is available for this script if you do not give arguments)
* Check for contamination with contamination-figure.R or the pandas-sample-similarity.ipynb notebook
    * Specify the files that were created by ConcatenateCloneFilesBatch.py
    * Specify the pt.table.csv that you got from the immunogenomics group
    * Check by hand if the column names in the pt.table are correct
    * Run the script
    * Usually I make reports for all samples per cell_type

## Preparation for Roche data (OBSOLETE) ##

* Use MakePTtableFromAAreads.R - Create a pt.table (sample description) from AA.reads file
* SplitAAreads.py - Splits the AA.reads table per sample (check the column names that you want to include in the file name!)
* SeqToFastq.py - Convert the tab-delimited files to fastq format

## How to cite

Barbera D. C. van Schaik, Paul L. Klarenbeek, Marieke E. Doorenspleet, Sabrina Pollastro, Anne Musters, Giulia Balzaretti, Rebecca E. Esveldt, Frank Baas, Niek de Vries and Antoine H. C. van Kampen (2016) T- and B-cell Receptor Repertoire Sequencing: Quality Control and Clone Identification. _In prep_.

## License

```
RESEDA - REpertoire SEquencing Data Analysis
Copyright (C) 2016-2020 Barbera DC van Schaik

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
```
