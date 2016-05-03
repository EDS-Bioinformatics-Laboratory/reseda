# WREP - Workflow for REPertoire sequencing

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

## How to run

The input files (in fastq format) can be specified by putting the paths in the file SAMPLES. At the top of execute-all.sh other parameters have to be set.

## How to cite

Barbera D. C. van Schaik, Paul L. Klarenbeek, Marieke E. Doorenspleet, Sabrina Pollastro, Anne Musters, Giulia Balzaretti, Rebecca E. Esveldt, Frank Baas, Niek de Vries and Antoine H. C. van Kampen (2016) T- and B-cell Receptor Repertoire Sequencing: Quality Control and Clone Identification. _In prep_.

## License
```
WREP - Workflow for REPertoire data analysis
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
