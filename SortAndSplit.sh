#!/bin/bash

myfiles=$@

for myfile in $myfiles; do
    sort ${myfile} > ${myfile}.sort
    split -l 14 ${myfile}.sort ${myfile}-
    rm ${myfile} ${myfile}.sort
done
echo "FINISHED"
