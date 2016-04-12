from __future__ import print_function
import sys
import os

if len(sys.argv)<2:
    sys.exit("Usage: create-weblogo.py all_info.csv (one or more files)")

myfiles = sys.argv[1:]

def getCDR3 (datafile):
    try:
        fh = open(datafile)
    except:
        sys.exit("cannot open file ", datafile)

    header = fh.readline().rstrip().split()
    c_cdr3 = header.index("cdr3pep")

    cdr3s = dict()
    for line in fh:
        line = line.rstrip()
        c = line.split()
        cdr3 = c[c_cdr3]
        cdr3_len = len(cdr3)
        cdr3s[cdr3_len] = cdr3s.get(cdr3_len, list())
        cdr3s[cdr3_len].append(cdr3)

    return(cdr3s)
        
def createSequenceLogo (outdir, cdr3_len, cdr3_list):
    inFASTA = outdir + "/seqlogo-cdr3-" + str(cdr3_len) + ".fasta"
    outPDF = outdir + "/seqlogo-cdr3-" + str(cdr3_len) + ".pdf"

    # Make CDR3 sequences unique
    cdr3_list = list(set(cdr3_list))
    fhOut = open(inFASTA, "w")
    for cdr3 in cdr3_list:
        print(cdr3, file=fhOut)
    fhOut.close()

    # Create the sequence logo
    mycmd = "./weblogo/seqlogo -F PDF -h 20 -w 40 -k 0 -c -S -f " + inFASTA + " > " + outPDF
    print(mycmd)
    os.system(mycmd)

################# MAIN ##############

for datafile in myfiles:
    outdir = "seqlogo-" + datafile.split("/")[-1]
    mycmd = "mkdir " + outdir
    print(mycmd)
    os.system(mycmd)

    cdr3s = getCDR3(datafile)
    n = 0
    for cdr3_len in sorted(cdr3s):
        # if n > 5:
        #     break
        createSequenceLogo(outdir, cdr3_len, cdr3s[cdr3_len]) # make sequence logo of CDR3
        n += 1
