from __future__ import print_function
import sys
from Bio import SeqIO

usage = "Usage: " + sys.argv[0] + " IMGT protein sequences with IMGT gaps"

if len(sys.argv) < 2:
    sys.exit(usage)

for myfile in sys.argv[1:]:
    try:
        fhIn = open(myfile, "r")
    except:
        sys.exit("cannot open file " + myfile)
    try:
        fhOut = open(myfile+".csv", "w")
    except:
        sys.exit("cannot create file " + myfile+".csv")

    print("V.gene,func,seq", file=fhOut)
    for record in SeqIO.parse(fhIn, "fasta") :
        c = record.description.split("|")
        (name,allele) = c[1].split("*")
        print(",".join([name,c[3],str(record.seq)]), file=fhOut)
