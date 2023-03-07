import pandas as pd
import gzip
import sys
from Bio import SeqIO

myfiles = sys.argv[1:]

for fastqFile in myfiles:
  fhOut = open(fastqFile.replace(".fastq.gz", "-ufis.csv"), "w")
  for record in SeqIO.parse(gzip.open(fastqFile, "rt"), "fastq"):
#  for record in SeqIO.parse(open(fastqFile), "fastq"):
    sequence = str(record.seq).upper()
    ufi = sequence[:6]
    print(record.id, ufi, sep="\t", file=fhOut)
  fhOut.close()

