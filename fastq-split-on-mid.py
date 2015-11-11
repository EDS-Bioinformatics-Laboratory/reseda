from __future__ import print_function
import sys
import os
import gzip
import regex
from sequences import *
from Bio import SeqIO

### Config ###

# indir = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/data"
# outdir = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/data/midsort"
indir = "."
outdir = "./midsort"

### Functions ###

def sortMIDS (motifs, fastq1File, fastq2File, outdir):
    '''
    In: motifs (list), fastq1File (file path), fastq2File (file path), outdir (directory path without slash at the end)
    Out: output files are generated in the directory specified in "outdir"
    '''

    prefixFastq1 = os.path.basename(fastq1File)
    prefixFastq1 = prefixFastq1.replace(".fastq.gz","")
    prefixFastq2 = os.path.basename(fastq2File)
    prefixFastq2 = prefixFastq2.replace(".fastq.gz","")

    print("Processing", prefixFastq1, "and", prefixFastq2)

    # Open fastq2, check the MIDs and create a barcode file
    try:
        fh2 = gzip.open(fastq2File, "rb")
    except:
        sys.exit("cannot open file:" + fastq2File)

    # Go through fastq entries of fastq2
    fh = dict()
    fh["nomatch"] = gzip.open(outdir + "/" + prefixFastq2 + "-nomatch.fastq.gz", "w")
    fh["report"] = open(outdir + "/" + prefixFastq2 + "-report.txt", "w")
    fh["midcount"] = open(outdir + "/" + prefixFastq2 + "-midcount.txt", "w")
    # i=0
    lookupTable = dict()
    midCount = dict()
    for record in SeqIO.parse(fh2, "fastq") :
        # if i>=10:
        #     break

        # Check which MID is present (4N MID 5nt)
        sequence = str(record.seq).upper()
        keepMatch = 0
        for motif in motifs:
            # Search for the motif
            m = regex.search(motif, sequence)
            if not m:
                next
            else:
                keepMatch = m

        
        # Write entry to the file with this particular MID
        # If there is no MID in the sequence write it to the file "nomatch"
        
        if keepMatch == 0:
            print(record.id, "nomatch", "nomatch", "nomatch", file=fh["report"])
            SeqIO.write(record, fh["nomatch"], "fastq")
            midCount["nomatch"] = midCount.get("nomatch", 0) + 1
        else:
            mid = keepMatch.group(2)
            if mid not in fh:  # If it is a new MID open a new output file
                fh[mid] = gzip.open(outdir + "/" + prefixFastq2 + "-" + mid + ".fastq.gz", "w")
            print(record.id, keepMatch.group(1), keepMatch.group(2), keepMatch.group(3), file=fh["report"])
            SeqIO.write(record, fh[mid], "fastq")
            lookupTable[record.id] = mid
            midCount[mid] = midCount.get(mid, 0) + 1

        # i = i + 1

    # Write mid counts to file
    for mid,freq in midCount.items():
        print(mid, freq, file=fh["midcount"])

    # Close all files
    for mykey in fh:
        fh[mykey].close()

    # Now open R1 and split that file based on the barcodes found in R2
    try:
        fh1 = gzip.open(fastq1File, "rb")
    except:
        sys.exit("cannot open file", fastq1File)

    fh = dict()
    fh["nomatch"] = gzip.open(outdir + "/" + prefixFastq1 + "-nomatch.fastq.gz", "w")
    # i=0
    for record in SeqIO.parse(fh1, "fastq") :
        # if i>=20:
        #     break

        if record.id not in lookupTable:
            SeqIO.write(record, fh["nomatch"], "fastq")
        else:
            mid = lookupTable[record.id]
            if mid not in fh:  # If it is a new MID open a new output file
                fh[mid] = gzip.open(outdir + "/" + prefixFastq1 + "-" + mid + ".fastq.gz", "w")
            SeqIO.write(record, fh[mid], "fastq")

        # i = i + 1

### Main ###

# Check if an argument was given to this script
# if len(sys.argv) < 2:
#     sys.exit('Usage: %s indir' % sys.argv[0])

# Input and configuration
midFile = "MIDS-miseq.txt"

# Read file with MIDs
motifs = readMotifsFromFile(midFile)

# Read 'indir' directory, store R1 and R2 files
fastq1File = list()
fastq2File = list()
for myfile in os.listdir(indir):
    if myfile.find("_R1_") > -1:
        fastq1File.append(myfile)
    elif myfile.find("_R2_") > -1:
        fastq2File.append(myfile)

# Create output directory if it doesn't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Sort the fastq files on MID
for i in range(len(fastq1File)):
    sortMIDS(motifs, indir+"/"+fastq1File[i], indir+"/"+fastq2File[i], outdir)

print("FINISHED")
