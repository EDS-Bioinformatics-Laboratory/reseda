Settings
========

Job monitoring
--------------

The samples can be divided over multiple (virtual) machines and be processed in
parallel. A lightweight job monitoring tool can be downloaded
here: https://bitbucket.org/barbera/progress

Settings webdav server
----------------------

These scripts assume that the data is on the SurfSara ResearchDrive (webdav server) and that the drive is mounted in /mnt/immunogenomics

* copy-from-beehub.sh
* copy-basespace-data-to-beehub.py
* execute-all.sh
* report-ALL.sh
* ConcatenateCloneFilesBatch.py
* MakeSamplesFiles.py
* VerifyBasespaceCopy.py

.. code-block:: bash

   sudo mount -t davfs -o uid=bioinfo,gid=bioinfo,rw https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics /mnt/immunogenomics

Create a .netrc file in your home directory:

.. code-block:: bash

   machine researchdrive.surfsara.nl
   login YOUR_RESEARCHDRIVE_USER
   password YOUR_PASSWORD
