#!/bin/bash

cdr3_files=`ls run*-vjcdr3-clones-*-IGH_HUMAN.csv`
match_scores="1 2 3 4 5"
mismatch_scores="-1 -2 -3 -4 -5"
gap_scores="-1 -2 -3 -4 -5"

for cdr3 in $cdr3_files; do
    for match in $match_scores; do
        for mismatch in $mismatch_scores; do
            for gap in $gap_scores; do
                echo "$cdr3 $match $mismatch $gap"
                python Dregion.py -cell IGH -match $match -mismatch $mismatch -gap $gap $cdr3
            done
        done
    done
done

