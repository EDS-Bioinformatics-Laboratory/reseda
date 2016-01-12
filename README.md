# Immunogenomics pipeline miseq

## Contributors

* Barbera van Schaik, b.d.vanschaik@amc.uva.nl
* Paul Klarenbeek, p.l.klarenbeek@amc.uva.nl
* Marieke Doorenspleet, m.e.doorenspleet@amc.uva.nl
* Sabrina Pollastro, s.pollastro@amc.uva.nl
* Guilia Balzaretti, g.balzaretti@amc.uva.nl
* Rebecca Esveldt, r.e.esveldt@amc.uva.nl
* Niek de Vries, n.devries@amc.uva.nl
* Antoine van Kampen, a.h.vankampen@amc.uva.nl

## Included third party software
* PEAR
* FastQC
* BWA
* VarScan
* Picard
* HLAforest

## To install
* Samtools
* Python 2.x, with libraries
    * __future__ import print_function
    * sys
    * os
    * subprocess
    * gzip
    * re
    * regex
    * biopython
    * editdistance
    * ?
* Perl
* Bioperl (for HLAforest)
* Bash

## How to run
See execute-all.sh

## Illumina data
With basemount data can be tranferred from the Illumina website
* http://blog.basespace.illumina.com/2015/07/23/basemount-a-linux-command-line-interface-for-basespace/
* mkdir basespace
* basemount basespace
