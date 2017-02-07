from __future__ import print_function
import sys
import numpy as np

try:
    fh0mm = open("/mnt/immunogenomics/RUNS/run06-20160306-miseq/results-tbcell/run06-first-amino-acid-0mm.txt")
    fh1mm = open("/mnt/immunogenomics/RUNS/run06-20160306-miseq/results-tbcell-1mm/run06-first-amino-acid-1mm.txt")
except:
    sys.exit("cannot open files")


def readFile(fh):
    readsList = dict()
    readsList["BCRh"] = list()
    readsList["TCRb"] = list()
    percList = dict()
    percList["BCRh"] = list()
    percList["TCRb"] = list()
    for line in fh:
        line = line.strip()
        (sample, aa, reads, perc) = line.split()

        if "BCRh" in sample:
            celltype = "BCRh"
        else:  # it is a TCRb sample
            celltype = "TCRb"

        if aa == "C":
            readsList[celltype].append(int(reads))
            percList[celltype].append(float(perc))

    for celltype in ["BCRh", "TCRb"]:
        print(celltype, sum(readsList[celltype]), np.mean(percList[celltype]))


readFile(fh0mm)
readFile(fh1mm)

fh0mm.close()
fh1mm.close()
