from __future__ import print_function
import sys
import os
import gzip
from Bio import SeqIO

def convertSffToFastq (sffFile):
    '''
    Description: convert sequences from sff to fasta format
    In: file.sff (str)
    Out: file.fastq.gz (str)
    '''
    fastqFile = os.path.basename(sffFile)
    fastqFile = fastqFile.replace(".sff",".fastq.gz")

    try:
        fhIn = open(sffFile)
        fhOut = gzip.open(fastqFile, "w")
    except:
        print("Could not open or write file:", sffFile, fastqFile)
        exit()

    for record in SeqIO.parse(fhIn, "sff"):
        SeqIO.write(record, fhOut, "fastq")


    fhIn.close()
    fhOut.close()

    return(fastqFile)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " sff-file(s)")

    myfiles = sys.argv[1:]
    for sffFile in myfiles:
        fastqFile = convertSffToFastq(sffFile)
        print("Converted:", sffFile, ">", fastqFile)
