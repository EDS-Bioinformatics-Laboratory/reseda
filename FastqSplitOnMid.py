from __future__ import print_function
import sys
import os
import gzip
import regex
from sequences import *
from Bio import SeqIO


def sortMIDS(motifs, fastqFile, outdir):
    '''
    In: motifs (list), fastqFile (file path), outdir (directory path without slash at the end)
    Out: output files are generated in the directory specified in "outdir"
    '''

    prefixFastq = os.path.basename(fastqFile)
    prefixFastq = prefixFastq.replace(".fastq.gz", "")

    print("Processing", prefixFastq)

    # Open fastq2, check the MIDs and create a barcode file
    try:
        fhIn = gzip.open(fastqFile, "rb")
    except:
        sys.exit("cannot open file:" + fastqFile)

    # Go through fastq entries of fastq file
    fh = dict()
    fh["nomatch"] = gzip.open(outdir + "/" + prefixFastq + "-nomatch.fastq.gz", "w")
    fh["report"] = open(outdir + "/" + prefixFastq + "-report.txt", "w")
    fh["midcount"] = open(outdir + "/" + prefixFastq + "-midcount.txt", "w")
    # i=0
    midCount = dict()
    for record in SeqIO.parse(fhIn, "fastq"):
        # if i>=10:
        #     break

        # Check which MID is present (4N MID 5nt)
        sequence = str(record.seq).upper()
        keepMatch = 0
        for motif in motifs:
            # Search for the motif in normal orientation
            m = regex.search(motif, sequence)
            if not m:
                # Search for the motif in reverse complement if it was not found
                m = regex.search(motif, comrev(sequence))
                if not m:
                    next
                else:
                    keepMatch = m
            else:
                keepMatch = m

        # Write entry to the file with this particular MID
        # If there is no MID in the sequence write it to the file "nomatch"

        if keepMatch == 0:
            print(record.id, "nomatch", "nomatch", "nomatch", file=fh["report"])
            SeqIO.write(record, fh["nomatch"], "fastq")
            midCount["nomatch"] = midCount.get("nomatch", 0) + 1
        else:
            # umi = keepMatch.group(1)
            # mid = keepMatch.group(2)
            # primerTB = keepMatch.group(3)
            umi = keepMatch.group(3)      # NEW
            mid = keepMatch.group(2)      # NEW
            primerTB = keepMatch.group(4) # NEW
            if mid not in fh:  # If it is a new MID open a new output file
                fh[mid] = gzip.open(outdir + "/" + prefixFastq + "-" + mid + ".fastq.gz", "w")
            print(record.id, umi, mid, primerTB, file=fh["report"])  # UMI, MID, Primer T or B
            SeqIO.write(record, fh[mid], "fastq")
            midCount[mid] = midCount.get(mid, 0) + 1

        # i = i + 1

    # Write mid counts to file
    for mid, freq in midCount.items():
        print(mid, freq, file=fh["midcount"])

    # Close all files
    for mykey in fh:
        fh[mykey].close()


if __name__ == '__main__':

    # Check if an argument was given to this script
    if len(sys.argv) < 4:
        sys.exit('Usage: %s midfile outdir fastq-file(s)' % sys.argv[0])

    [midFile, outdir] = sys.argv[1:3]
    fastqFiles = sys.argv[3:]

    # midFile = "MIDS-miseq.txt"
    # indir = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/data"
    # outdir = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/data/midsort"

    # Read file with MIDs
    motifs = readMotifsFromFile(midFile)

    # Create output directory if it doesn't exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Sort the fastq files on MID
    for i in range(len(fastqFiles)):
        sortMIDS(motifs, fastqFiles[i], outdir)

    print("FINISHED")
