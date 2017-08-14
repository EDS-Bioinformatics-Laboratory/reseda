# RESEDA - REpertoire SEquencing Data Analysis

Data analysis workflow for T- and B-cell receptor repertoire sequencing.
The workflow identifies clones and their frequency from next generation sequencing of repertoires and includes steps for quality control and bias correction.

## Workflow

![workflow](workflow.png)

## Required software

* PEAR
* FastQC
* BWA
* VarScan
* Picard
* Samtools

## Other requirements

* Bash
* Python 2.7
    * future (print_function)
    * sys
    * os
    * subprocess
    * gzip
    * re
    * regex
    * math
    * random
    * sqlite3
    * matplotlib
    * biopython
    * editdistance
    * alignment
    * numpy
    * scipy
* R
    * beeswarm

## Job monitoring

Tip: divide the samples over multiple (virtual) machines and run everything in parallel. You can download a lightweight job monitoring tool [HERE](https://bitbucket.org/barbera/progress).

## How to run - Using only the code in this repository

The input files (in fastq format) can be specified by putting the paths in the file SAMPLES. At the top of execute-all.sh you see which other parameters have to be set.

Example: ./execute-all.sh output-dir-on-webdav-server MIDS-miseq.txt human IGH IGH_HUMAN

## How to run - Using Topos for sending jobs and a job monitor tool

### Preparation ###
* Mount basespace. Instructions are in basespace.txt
* Specify the basespace directories in copy-basespace-data-to-beehub.py and run: python copy-basespace-data-to-beehub.py > basespace-copy-data.sh
* Verify if all data has been copied with VerifyBasespaceCopy.py
* Convert the MiSeq sample sheet with MetaData.py (creates a json file)
* Mount the beehub webdav server
* Add extra information to the json file with MakeSamplesFiles.py (this will also make the SAMPLE-* files)
* Make manageable jobs by splitting the SAMPLE-* files, e.g.: split -l 20 SAMPLES-run13-human-BCRh SAMPLES-run13-human-BCRh-
* Check manually if the nr of lines in the SAMPLE-* files (the total divided by 2) is equal to the nr of samples in the sample sheet
* Create Topos jobs with ToposCreateTokens.py
* Upload Topos jobs with ToposUploadFiles.py

### Starting the jobs ###
* Start virtual machines for the analysis (in the SurfSara HPC cloud webinterface)
* Add all ip-adresses to the ip-list file in the '../progress' directory
* Transfer the setup-and-run.sh by running roll-out-scripts.py and executing the code that was generated by this roll-out script
* Login to each machine and run ./setup-and-run.sh

### When all jobs are finished ###
* In the jobs the data is automatically transferred to the beehub webdav server
* Specify the mounted directory and sample.json file in ConcatenateCloneFilesBatch.py and run it to concatenate the clone files per project+organism+cell_type
* Specify the mounted directory in report-ALL.sh and run it to generate reports about the sequence run
* Check for contamination with contamination-figure.R
    * Specify the files that were created by ConcatenateCloneFilesBatch.py
    * Specify the pt.table.csv that you got from the immunogenomics group
    * Check by hand if the column names in the pt.table are correct
    * Run the script
    * Usually I make reports for all samples per project+cell_type and one report for all the samples in a run

## How to cite

Barbera D. C. van Schaik, Paul L. Klarenbeek, Marieke E. Doorenspleet, Sabrina Pollastro, Anne Musters, Giulia Balzaretti, Rebecca E. Esveldt, Frank Baas, Niek de Vries and Antoine H. C. van Kampen (2016) T- and B-cell Receptor Repertoire Sequencing: Quality Control and Clone Identification. _In prep_.

## License
```
RESEDA - REpertoire SEquencing Data Analysis
Copyright (C) 2016 Barbera DC van Schaik

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