#!/bin/bash

POOLNAME=d8c24f78f9772cbdff54cf62

# Mount beehub
./mount-beehub.sh
wait

# Start progess monitor
cd git/progress
wait
screen ./run.sh
wait

# Get latest version of the pipeline
cd ../tbcell-miseq-pipeline/
wait
git pull origin devel
wait

# Put reference files in the right place
mv reference/* .
wait

# Start analysis
nohup python ToposGetTokens.py ${POOLNAME} > nohup.out 2> nohup.err < /dev/null &
