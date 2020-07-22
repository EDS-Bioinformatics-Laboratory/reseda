Run with PiCaS
==============

How to run - Using PiCaS for sending jobs and using a job monitor tool

Note: this is what Barbera does for each sequence run

PiCaS: https://picas.surfsara.nl:6984/_utils/database.html?tbcellrep-bschaik-hpc

Job monitoring tool: https://bitbucket.org/barbera/progress/

Helper scripts
--------------

Fill in the run name and output directory in ENV.sh and source the script. You
can use the environment variables when you go through the steps below.

.. code-block:: bash

   RUN=run30-20180723-miseq
   OUTDIR=results-tbcell
   WEBDAV=https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/${RUN}
   POOLNAME=d8c24f78f9772cbdff54cf62

Files can be copied from the webdav server by adding file paths in SAMPLES and
then use:

.. code-block:: bash

   ./copy-from-beehub.sh

Example of the content in the SAMPLES file:

.. code-block:: bash

   /mnt/immunogenomics/RUNS/TESTDATA/data/test-B_S1_L001_R1_001.fastq.gz
   /mnt/immunogenomics/RUNS/TESTDATA/data/test-B_S1_L001_R2_001.fastq.gz
   /mnt/immunogenomics/RUNS/TESTDATA/data/test-T_S2_L001_R1_001.fastq.gz
   /mnt/immunogenomics/RUNS/TESTDATA/data/test-T_S2_L001_R2_001.fastq.gz

Data can be transferred from the webdav server with:

.. code-block:: bash

   ./copy-to-webdav.sh WEBDAVURL FILE1 [FILE2]

Preparation
-----------

Mount the ResearchDrive from the home directory on the cloud machine and create a directory from the new run.

.. code-block:: bash

    ./mount-beehub.sh
    mkdir /mnt/immunogenomics/RUNS/run35-20190609-miseq
    mkdir /mnt/immunogenomics/RUNS/run35-20190609-miseq/data

Mount basespace from the home directory. Instructions are in basespace.txt (you need the password for Illumina base space, Rebecca has it) The latest sequence run is in basespace/Projects/

.. code-block:: bash

    basemount basespace

Specify the run and the basespace sub-directories as argument to
copy-basespace-data-to-beehub.py and run it. The file basespace-copy-data.sh
and basespace-calc-checksum.sh will be created. Run these scripts to copy the
data from basespace to the ResearchDrive and to calculate the SHA1 sums
(last one is needed for the VerifyBasespaceCopy.py script)

.. code-block:: bash

   python copy-basespace-data-to-beehub.py -r runNN-yyyymmdd-miseq run-dir-in-project-dir-on-basespace
   nohup bash basespace-copy-data.sh > nohup-copy.out 2> nohup-copy.err < /dev/null &
   nohup bash basespace-calc-checksum.sh > nohup-check.out 2> nohup-check.err < /dev/null &

The immunogenomics group normally sends a semicolon-separated files with all the required information. When you receive a so-called 'pt-table' file you need to copy/paste this information in a 'datasheet'. Download one of the earlier datasheets (e.g. from the previous run), edit the header to match it with the latest sequence run and copy/paste the right columns in this sheet. The order of the columns is not important, but the column names are.

Convert the MiSeq sample sheet (datasheet) with MetaData.py (creates a json file)

.. code-block:: bash

   python MetaData.py Miseq-sample-Datasheet.csv > Miseq-sample-Datasheet.json

Mount the ResearchDrive webdav server if you have not done so already.

.. code-block:: bash

   sudo mount -t davfs -o uid=bioinfo,gid=bioinfo,rw https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics /mnt/immunogenomics

Verify if all data has been copied (optional). Note that the SHAsum verification does not work currently.

.. code-block:: bash

   python2 VerifyBasespaceCopy.py -i yyyymmdd_RUNxx_Datasheet.json -r runNN-yyyymmdd-miseq

You need the mounted ResearchDrive directory name and the json file from the previous step.

Add extra information to the json file with MakeSamplesFiles.py (this will also make the SAMPLE-* files)

.. code-block:: bash

   python MakeSamplesFiles.py -r yyyymmdd-RUNnn-datasheet.json -w /mnt/immunogenomics/RUNS/runNN-yyyymmdd-miseq/data/

Sort and split the SAMPLE-* files with: SortAndSplit.sh It does the following:

    * Sorts the SAMPLE-* files: sort SAMPLE-blah > SAMPLE-blah.sort
    * Makes manageable jobs by splitting the SAMPLE-\*.sort files, e.g.: split -l 20 SAMPLES-run13-human-BCRh.sort SAMPLES-run13-human-BCRh-

.. code-block:: bash

   ./SortAndSplit.sh SAMPLE-*

Create PiCaS jobs with ToposCreateTokens.py (run with the -h option to see the arguments)

.. code-block:: bash

   python ToposCreateTokens.py -r runNN-yyyymmdd-miseq -m MIDS-miseq-umi.txt -o results-tbcell -p paired -b yes -u yes

Upload PiCaS jobs with picas/createTokens.py JSON-FILES (json files that were created in the previous step)

.. code-block:: bash

   cd picas
   python createTokens.py ../tokens/*

Starting the jobs
-----------------

Start virtual machines for the analysis (in the SurfSara HPC cloud webinterface)
https://ui.hpccloud.surfsara.nl/

Add all ip-adresses to the ip-list file in the '../progress' directory (the job monitoring tool)

Transfer the setup-and-run.sh by running roll-out-scripts.py and executing the
code that was generated by this roll-out script

.. code-block:: bash

   python roll-out-scripts.py > tmp.sh
   bash tmp.sh

Login to each machine and run:

.. code-block:: bash

   ./setup-and-run.sh

When you see "Serving HTTP on 0.0.0.0 port 8000 ..." do ctrl+a d. This will detach the current "screen" and run the process in the background. After this the script will continue. Log out when you get the prompt back. The ToPoS job has started in the background. Repeat for all the other virtual machines.

Check progress of the analysis
------------------------------

Go to the "progress" repository and run the display.py script. When all jobs are "FINISHED" you can move on and generate clone files and reports.

.. code-block:: bash

    cd ../progress
    ./display.py

The status of the jobs can also be traced via PiCas. URL-HERE The analysis is finished when there are no more jobs (tokens) left.

Monitoring PiCaS via the command-line:

.. code-block:: bash

    cd picas
    python createViews.py

When all jobs are finished
--------------------------

In the jobs the data is automatically transferred to the ResearchDrive webdav server

Execute ConcatenateCloneFilesBatch.py to generate a bash script for
concatenating clone files per project+organism+cell_type (run the generated script)

.. code-block:: bash

   python ConcatenateCloneFilesBatch.py -r yyyymmdd-RUNnn-datasheet-new.json -w /mnt/immunogenomics/RUNS/runNN-yyyymmdd-miseq/results-tbcell/final/ > tmp.sh
   nohup bash tmp.sh > nohup-concat.out 2> nohup-concat.err < /dev/null &


Run report-ALL.sh to generate reports about the sequence run (help is available for this script if you do not give arguments)

.. code-block:: bash

   ./report-ALL.sh -r runNN-YYYYMMDD-miseq -i yyyymmdd_RUNnn_Datasheet-new.json -b yes -o results-tbcell

Check for contamination with contamination-figure.R or the pandas-sample-similarity.ipynb notebook

* Specify the files that were created by ConcatenateCloneFilesBatch.py
* Specify the pt.table.csv that you got from the immunogenomics group
* Check by hand if the column names in the pt.table are correct
* Run the script
* Usually the reports are made for all samples per cell_type

Starting the jupyter notebook server:

.. code-block:: bash

   jupyter notebook --no-browser

The (password-protected) notebook will be accessible via the ip address of the virtual machine and port 8888. E.g. http://145.100.57.10:8888/
