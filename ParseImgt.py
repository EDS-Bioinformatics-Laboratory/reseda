from __future__ import print_function
from Bio import SeqIO

# fh = open("imgt-record.txt")
fh = open("/home/barbera/Downloads/imgt.dat")
for record in SeqIO.parse(fh, "imgt"):
    for f in record.features:
        # print(f)
        # print("--------")
        print("LOCATION IS:", f.location)
        print("--------")

fh.close()
