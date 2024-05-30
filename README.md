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

* Python 3.5 (or higher)
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

### Settings for curl

Create a .netrc file in your home directory:

* machine researchdrive.surfsara.nl
* login YOUR_RESEARCHDRIVE_USER
* password YOUR_PASSWORD

### Mount ResearchDrive
```
# change "bioinfo" with your own user name on the machine you are working on
sudo mount -t davfs -o uid=bioinfo,gid=bioinfo,rw https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics /mnt/immunogenomics
```

### Rclone installation and configuration

You can download the deb package, or fetch it from the [ansible-playbooks](https://github.com/EDS-Bioinformatics-Laboratory/Ansible-Playbooks/tree/master/Reseda) git repository

``sudo apt install ./rclone-v1.61.1-linux-amd64.deb``

Run: ``rclone config``

```
name> remote
Storage> webdav
url> https://researchdrive.surfsara.nl/remote.php/nonshib-webdav
vendor> owncloud
```

## How to run the code

You have received the raw data (fastq files) and sample information from the Immunogenomics group.
These are stored on the ResearchDrive. Follow the steps below to store it on the right location.

You might need to convert the sample information into the Miseq Datasheet format first.
This can be done with the notebook ``MakeDatasheetFromPT.ipynb``.

### Setting up the FSS data structure and create a new git branch

Use the script from the ENCORE_AUTOMATION repository:

https://github.com/EDS-Bioinformatics-Laboratory/ENCORE_AUTOMATION/blob/main/Processing/2-CREATE-TEMPLATE-REPSEQ/README.md

Transfer the fastq files to the appropriate directory on the ResearchDrive:
```
/mnt/immunogenomics/RUNS/runNN-yyyymmdd-miseq/Data/NameOfDataset_1/Raw/
```

### Preparation ###

* Start virtual machines for the analysis (in the [SurfSara cloud webinterface](https://portal.live.surfresearchcloud.nl/))
* An Ansible script is used to install the machines
    * See: [../ansible-playbooks/Reseda/README.md](../ansible-playbooks/Reseda/README.md)

### Create tokens (jobs) and upload them to the virtual machines

Note: the following steps can run on a Linux laptop or on a virtual machine. In the latter case you need to install the software on the VM first (with Ansible)

* Convert the MiSeq sample sheet with ``MetaData.py`` (creates a json file)
    * An example of a datasheet can be found in the run directories on the ResearchDrive (Data/NameOfDataset_1/Meta/)
* Mount the ResearchDrive webdav server, if you haven't done so already
* Add extra information to the json file with ``MakeSamplesFiles.py`` (this will also make the SAMPLE-* files)
* Sort and split the SAMPLE-* files with: ``./SortAndSplit.sh SAMPLE-*`` It does the following:
    * Sorts the SAMPLE-* files: sort SAMPLE-blah > SAMPLE-blah.sort
    * Makes manageable jobs by splitting the SAMPLE-\*.sort files, e.g.: split -l 20 SAMPLES-run13-human-BCRh.sort SAMPLES-run13-human-BCRh-
* Create jobs with ``ToposCreateTokens.py`` (run with the ``-h`` option to see the arguments)
* Upload the jobs to the VMs using the script ``TransferTokens.py`` (run with ``-h`` to see the options)

### Start the jobs on each virtual machine

* Login to each virtual machine and do the following:
  * Copy/move the relevant (or all) database files from directories "reference", "reftables" and "mids" to the root directory of the repository: ``git/reseda/``
  * Start the script ``RUN-RESEDA.py tokens/``

### When all jobs are finished ###

* In the jobs the results are automatically transferred to the ResearchDrive webdav server
* Check with the notebook ``CompareProcessedUnprocessed.ipynb`` if all result files are on the ResearchDrive
* Execute ``ConcatenateCloneFilesBatch.py`` to generate a bash script for concatenating clone files per project+organism+cell_type (run the generated script)
* Run ``report-ALL.sh`` to generate reports about the sequence run (help is available for this script if you do not provide arguments to the script)
* Check for contamination with the notebooks ``SampleSimilarity.ipynb`` and ``SharedClonesDirection.ipynb``
    * Specify the files that were created by ConcatenateCloneFilesBatch.py
    * Specify the Datasheet table that you received from the immunogenomics group
    * Check by hand if the column names in the Datasheet are correct
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
Copyright (C) 2016-2024 Barbera DC van Schaik

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
