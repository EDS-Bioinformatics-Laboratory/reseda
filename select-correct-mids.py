from __future__ import print_function
import sys

# Create this file with: wc -l *-all_info.csv > wc-run04.txt
wcFile = sys.argv[1]   # "wc.txt"

try:
    fh = open(wcFile, "r")
except:
    sys.exit("cannot open wcFile")

# Store wc in this dictionary structure: d[(sample, celltype)][mid] = wc

d = dict()
for line in fh:
    line = line.strip()  # remove whitespaces before and after line
    (wc, filepath) = line.split() # get wordcount and filepath
    wc = int(wc)
    path = filepath.split("/")    # get directories and filename

    filename = path[-1]
    elements = filename.split("-") # filename contains info about:
    celltype = elements[-2]
    mid = elements[-3]
    sample = "-".join(elements[:-3])

    if mid != "nomatch":
        d[(sample,celltype)] = d.get((sample,celltype), dict()) # create new dict if not exist
        d[(sample,celltype)][mid] = wc

for (sample,celltype) in sorted(d):
    d_tmp = d[(sample, celltype)]
    topmid = sorted(d_tmp, key=d_tmp.get, reverse=True)[0]
    #print("Highest MID:", sample, celltype, topmid, d[(sample,celltype)][topmid])
    print("mv", sample+"-"+topmid+"*", "correct-mid/")
