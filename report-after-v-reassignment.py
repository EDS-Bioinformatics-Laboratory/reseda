from __future__ import print_function
import sys

# Input
if len(sys.argv) < 2:
    sys.exit("Usage: report-after-v-reassignment.py /path/to/*.rr.clones_subs.csv")
datafiles = sys.argv[1:]


def countReads(datafile):
    try:
        fh = open(datafile)
    except:
        sys.exit("cannot open file " + datafile)

    fh.readline()
    totalreads = 0
    for line in fh:
        line = line.strip()
        c = line.split()
        freq = int(c[3])
        totalreads += freq

    fh.close()

    return(totalreads)


try:
    fhOut = open("report-AFTER-V-REASSIGNMENT.txt", "w")
except:
    sys.exit("cannot write to file")

for datafile in datafiles:
    totalreads = countReads(datafile)
    print(datafile, totalreads)
    print(datafile, totalreads, file=fhOut)

fhOut.close()
