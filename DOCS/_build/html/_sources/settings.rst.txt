Settings
========

Job monitoring
--------------

The samples can be divided over multiple (virtual) machines and be processed in
parallel. A lightweight job monitoring tool can be downloaded
here: https://bitbucket.org/barbera/progress

Usernames and passwords
-----------------------

The following accounts are necessary:

* Illumina basespace (via Rebecca Esveldt or Barbera van Schaik)
* HPC cloud account with access to the 'tbcellrep-amc' group
* Barbera needs to add your public SSH key to the virtual image (needs to be done once)
* 'bioinfo' password on the cloud machine
* Your own ResearchDrive account with access to the immunogenomics repository
* Password of the jupyter notebook server

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

Settings cloud environment
--------------------------

For security reasons one can only login to the cloud virtual machines with public-private key pairs.
This can normally be configured in the HPC cloud dashboard https://ui.hpccloud.surfsara.nl/
via your-user-name > Settings > Public SSH key

For the immunogenomics machines this is unfortunately not the case.
Because the machines have been installed a while ago and everything runs under the user 'bioinfo' your public key needs to be added by Barbera once.
She will add your public key to the .ssh/authorized_keys file on the cloud machine and save this as a new image.
