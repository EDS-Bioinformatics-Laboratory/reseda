from __future__ import print_function
import sys
import os
import gzip
from Bio import SeqIO


def convertSffToFastq(f):
    '''
    Description: convert sequences from sff to fasta format
    In: file.sff (str)
    Out: file.fastq.gz (str)
    '''
    fastqFile = os.path.basename(f)
    fastqFile = fastqFile.replace(".sff", ".fastq.gz")

    try:
        fhIn = open(f)
        fhOut = gzip.open(fastqFile, "w")
    except:
        print("Could not open or write file:", f, fastqFile)
        exit()

    for record in SeqIO.parse(fhIn, "sff"):
        SeqIO.write(record, fhOut, "fastq")

    fhIn.close()
    fhOut.close()

    return(fastqFile)


def convertFastaToFastq(f):
    '''
    Description: convert sequences from sff to fasta format
    In: file.sff (str)
    Out: file.fastq.gz (str)
    '''
    fastqFile = os.path.basename(f)
    fastqFile = fastqFile.replace(".fasta", ".fastq.gz")

    try:
        fhIn = open(f)
        fhOut = gzip.open(fastqFile, "w")
    except:
        print("Could not open or write file:", f, fastqFile)
        exit()

    for record in SeqIO.parse(fhIn, "fasta"):
        record.letter_annotations["phred_quality"] = [40] * len(str(record.seq))
        SeqIO.write(record, fhOut, "fastq")

    fhIn.close()
    fhOut.close()

    return(fastqFile)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " sequence-file(s) Supported formats: sff and fasta")

    myfiles = sys.argv[1:]
    for f in myfiles:
        fastqFile = convertFastaToFastq(f)
        print("Converted:", f, ">", fastqFile)
