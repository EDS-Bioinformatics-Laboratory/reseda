#!/bin/bash

# Put reference files in the work directory
mv reference/* .
wait

# Put MIDS files in the work directory
mv mids/* .
wait

# Put the ref.table.* files in the work directory
mv reftables/* .
wait

nohup python RUN-RESEDA.py tokens/ > nohup.out 2> nohup.err < /dev/null &

