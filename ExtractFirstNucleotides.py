from __future__ import print_function
import sys
import gzip
from Bio import SeqIO


def getFirstNucleotides(myfile):
    fh = gzip.open(myfile, "rt")
    for record in SeqIO.parse(fh, "fastq"):
        seq = str(record.seq)
        print(record.id, seq[:5], seq[:6], seq[:7], seq[:8], seq[:9], seq[:10], seq[:11], seq[:12], seq[:13], seq[:14], sep="\t")
    fh.close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:", sys.argv[0], "fastq-file(s)")
        exit()

    print("acc\tnt5\tnt6\tnt7\tnt8\tnt9\tnt10")

    for myfile in sys.argv[1:]:
        getFirstNucleotides(myfile)
