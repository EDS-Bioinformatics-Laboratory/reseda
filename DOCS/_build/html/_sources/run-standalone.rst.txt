Run standalone
==============

How to run, using only the code in this repository

Copy/move the relevant (or all) database files from directory "reference" to the
root directory of the repository (same directory as the execute-all.sh script)

The input files (in fastq format) can be specified by putting the paths in the
file SAMPLES. Running ./execute-all.sh without arguments shows you which
parameters can be set.

Example:

.. code-block:: bash

   ls TESTDATA/test-B*.fastq.gz > SAMPLES
   ./execute-all.sh -r mytestrun -l local -m MIDS-miseq.txt -org human -cell IGH -celltype IGH_HUMAN -u no

The results will be on the machine where you run this script, and if the webdav
settings are configured well a new directory with all the files will appear on
the ResearchDrive:
https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/mytestrun

MIDs: You need MIDS-miseq.txt for samples without UMI and MIDS-miseq-umi.txt for samples with UMI
