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


def convertFastaQualToFastq(f,q):
    '''
    Description: convert fasta+qual files to fastq
    In: file.fasta and file.qual (both str)
    Out: file.fastq.gz (str)
    '''

    fastqFile = os.path.basename(f)
    fastqFile = fastqFile.replace(".seq", ".fastq.gz")
    try:
        fhOut = gzip.open(fastqFile, "w")
    except:
        print("Could not write file:", fastqFile)
    reads = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    for rec in SeqIO.parse(q, "qual"):
        reads[rec.id].letter_annotations["phred_quality"]=rec.letter_annotations["phred_quality"]
        SeqIO.write(reads[rec.id], fhOut, "fastq")
    fhOut.close()

    return(fastqFile)


def convertFastqToFasta(f):
    '''
    Description: convert sequences from fastq to fasta format
    In: file.fastq.gz (str)
    Out: file.fasta (str)
    '''
    fastaFile = os.path.basename(f)
    fastaFile = fastaFile.replace(".fastq.gz", ".fasta")

    try:
        fhIn = gzip.open(f)
        fhOut = open(fastaFile, "w")
    except:
        print("Could not open or write file:", f, fastaFile)
        exit()

    for record in SeqIO.parse(fhIn, "fastq"):
        SeqIO.write(record, fhOut, "fasta")

    fhIn.close()
    fhOut.close()

    return(fastaFile)


def convertTabToFastq(f):
    '''
    Description: convert sequences from tab to fasta format
    In: file.txt (str)
    Out: file.fastq.gz (str)
    '''
    fastqFile = os.path.basename(f)
    fastqFile = fastqFile.replace(".txt", ".fastq.gz")

    try:
        fhIn = open(f)
        fhOut = gzip.open(fastqFile, "w")
    except:
        print("Could not open or write file:", f, fastqFile)
        exit()

    for record in SeqIO.parse(fhIn, "tab"):
        record.letter_annotations["phred_quality"] = [40] * len(str(record.seq))
        SeqIO.write(record, fhOut, "fastq")

    fhIn.close()
    fhOut.close()

    return(fastqFile)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit("Usage: " + sys.argv[0] + " sff|fasta|fastaqual|tab|fastq2fasta sequence-file(s) Supported formats: sff, fasta, fastaqual, tab, fastq2fasta")

    input_type = sys.argv[1]
    myfiles = sys.argv[2:]
    for i, f in enumerate(myfiles):
        if input_type == "fasta":
            fastqFile = convertFastaToFastq(f)
        elif input_type == "sff":
            fastqFile = convertSffToFastq(f)
        elif input_type == "fastq2fasta":
            fastqFile = convertFastqToFasta(f)
        elif input_type == "tab":
            id_col = 0
            seq_col = 1
            fastqFile = convertTabToFastq(f)
        elif input_type == "fastaqual":
            q = myfiles[i + 1]
            fastqFile = convertFastaQualToFastq(f,q)
            print("Converted:", f, q, ">", fastqFile)
            exit(0)
        else:
            sys.exit("Unknown input type: " + input_type)
        print("Converted:", f, ">", fastqFile)
