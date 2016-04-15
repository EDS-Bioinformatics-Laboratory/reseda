from __future__ import print_function
import sys

if len(sys.argv) < 2:
    sys.exit("Usage: concatenate-clone-files.py *clones_subs.csv")

myfiles = sys.argv[1:]

# Open file for writing
try:
    fhOut = open("run-clones_subs.csv", "w")
except:
    sys.exit("cannot write to file")

########### FUNCTIONS ##########

def parseSampleName (myfile):
    '''
    Description: extract the sample name from the file name
    In: file name (str)
    Out: sample name (str)
    '''
    path = myfile.split("/")[-1]
    sample, rest = path.split("_L001.assembled-")
    mid = rest.split("-")[0]
    return(sample,mid)

def printContent (myfile, fhOut):
    '''
    Description: open file, skip header and print to file
    '''
    try:
        fh = open(myfile)
    except:
        sys.exit("cannot open file "+myfile)

    sample,mid = parseSampleName(myfile)

    header = fh.readline() # skip header
    for line in fh:
        line = [sample,mid] + line.strip().split()
        print("\t".join(line), file=fhOut)

    fh.close()

############ MAIN ############

# Read header of the first file and write to disk
try:
    fh = open(myfiles[0])
except:
    sys.exit("cannot open file "+myfiles[0])

header = fh.readline()
header = ["Sample","MID"] + header.strip().split()
print("\t".join(header), file=fhOut)
fh.close()

# Read content of all files and write to fhOut
for myfile in myfiles:
    printContent(myfile, fhOut)

fhOut.close()

print("DONE")
