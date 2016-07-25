from __future__ import print_function
import sys
import json

# Usage: python ToposCreateTokens.py runNN-YYYYMMDD-miseq MIDS-miseq.txt SAMPLE-files

# Create pool: topos newPool
# Upload tokens: uploadFileAsToken [POOLNAME] [FILENAME]                                  #
#            or: uploadFilesInDirAsTokens [POOLNAME] [DIRNAME] 
# Token pool is here: https://topos.grid.sara.nl/4.1/pools/POOLNAME/

def getSamples (myfile):
    '''
    Description: read sample names and return a list of all the file names
    In: SAMPLE-file
    Out: list
    '''
    samples = list()
    try:
        fh = open(myfile)
    except:
        sys.exit("cannot open file: " + myfile)

    for line in fh:
        line = line.strip()
        samples.append(line)

    fh.close()
    return(samples)

def writeJson (f, js):
    '''
    Description: write json code to a file
    In: f (filename str), js (json dict)
    Out: -
    '''
    try:
        fhOut = open(f, "w")
    except:
        sys.exit("cannot create file: " + f)

    print(json.dumps(js, indent=4), file=fhOut)

    fhOut.close()
    print("Wrote", f, "to disk")

def guessCellAndOrganism (myfile):
    '''
    Description: the file name usually contains the celltype and organism, guess which one it is
    In: filename (str)
    Out: cell (str), organism (str), celltype (str)
    '''
    cell = "UNKNOWN"
    organism = "UNKNOWN"

    if "human" in myfile:
        organism = "human"
    elif "mouse" in myfile:
        organism = "mouse"

    if "IGH" in myfile or "BCRh" in myfile:
        cell = "IGH"
    elif "TCRb" in myfile:
        cell = "TRB"
    elif "TCRa" in myfile:
        cell = "TRA"
    elif "IGL" in myfile:
        cell = "IGL"
    elif "IGK" in myfile:
        cell = "IGK"
    elif "HLA" in myfile:
        cell = "HLA"

    celltype = cell + "_" + organism.upper()
    return(cell, organism, celltype)

if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("Usage: " + sys.argv[0] + " runNN-YYYYMMDD-miseq MIDS-miseq.txt SAMPLE-files")

    run = sys.argv[1]
    mids = sys.argv[2]

    for myfile in sys.argv[3:]:
        outfile = "tokens/" + myfile.split("/")[-1] + ".json"
        (cell, organism, celltype) = guessCellAndOrganism(myfile)
        js = {
            "run": run,
            "cell": cell,
            "organism": organism,
            "celltype": celltype,
            "mids": mids
        }
        js["samples"] = getSamples(myfile)
        writeJson(outfile, js)
