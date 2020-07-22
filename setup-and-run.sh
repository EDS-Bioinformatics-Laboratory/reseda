#!/bin/bash

POOLNAME=d8c24f78f9772cbdff54cf62

# Mount beehub
./mount-beehub.sh
wait

# Clean up the apt cache
sudo apt clean

# Start progess monitor
cd git/progress
wait
git pull origin master
wait
screen ./run.sh
wait

# Get latest version of the pipeline
cd ../reseda/
wait
git pull origin master
wait

# Put reference files in the work directory
mv reference/* .
wait

# Put MIDS files in the work directory
mv mids/* .
wait

# Put the ref.table.* files in the work directory
mv reftables/* .
wait

# Put all scripts necessary for Picas in the work directory
mv picas/* .
tar -xvf picas.tar
tar -zxf couchdb.tgz

# This is temporary until the library is installed on the virtual machine
pip install couchdb

# Start analysis
#nohup python ToposGetTokens.py ${POOLNAME} > nohup.out 2> nohup.err < /dev/null &
echo "Start the pilot job tasks by contacting PiCaS tokens"
python2 pilot.py
