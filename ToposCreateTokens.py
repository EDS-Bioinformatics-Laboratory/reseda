from __future__ import print_function
import sys
import json
import argparse

# Usage: python ToposCreateTokens.py runNN-YYYYMMDD-miseq MIDS-miseq.txt SAMPLE-files

# Create pool: topos newPool
# Upload tokens: uploadFileAsToken [POOLNAME] [FILENAME]                                  #
#            or: uploadFilesInDirAsTokens [POOLNAME] [DIRNAME]
# Token pool is here: https://topos.grid.sara.nl/4.1/pools/POOLNAME/


def getSamples(myfile):
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


def writeJson(f, js):
    '''
    Description: write json code to a file
    In: f (filename str), js (json dict)
    Out: -
    '''
    try:
        fhOut = open(f, "w")
    except:
        sys.exit("cannot create file: " + f)

    print(json.dumps(js), file=fhOut)

    fhOut.close()
    print("Wrote", f, "to disk")


def guessCellAndOrganism(myfile):
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

    if "IGH" in myfile or "BCRh" in myfile or "IgM" in myfile or "IgG" in myfile or "IgD" in myfile:
        cell = "IGH"
    elif "TCRb" in myfile:
        cell = "TRB"
    elif "TCRa" in myfile:
        cell = "TRA"
    elif "IGL" in myfile or "BCRl" in myfile:
        cell = "IGL"
    elif "IGK" in myfile or "BCRk" in myfile:
        cell = "IGK"
    elif "HLA" in myfile:
        cell = "HLA"

    celltype = cell + "_" + organism.upper()
    return(cell, organism, celltype)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create ToPoS tokens (json)')
    parser.add_argument('-r', '--run', default='runNN-YYYYMMDD-miseq', type=str, help='Sequence run directory (default: %(default)s)')
    parser.add_argument('-m', '--mids', default='MIDS-miseq-umi.txt', type=str, help='MID scheme file name (default: %(default)s)')
    parser.add_argument('-o', '--outdir', default='results-tbcell', type=str, help='Output directory (default: %(default)s)')
    parser.add_argument('-p', '--protocol', default='paired', type=str, help='single or paired (default: %(default)s)')
    parser.add_argument('-b', '--barcodes', default="yes", type=str, help='Were additional internal MIDs used? yes/no (default: %(default)s)')
    parser.add_argument('-u', '--umis', default="yes", type=str, help='Were UMIs used? yes/no (default: %(default)s)')
    parser.add_argument('-mm', '--mismatches', default=0, type=str, help='Number of mismatches in CDR3 extraction (default: %(default)s)')
    parser.add_argument("sample_files", type=str, nargs='+', help='Path(s) to SAMPLE file(s)')
    args = parser.parse_args()

    if args.run == "runNN-YYYYMMDD-miseq":
        parser.print_help()
        exit()

    for myfile in args.sample_files:
        outfile = "tokens/" + myfile.split("/")[-1] + ".json"
        (cell, organism, celltype) = guessCellAndOrganism(myfile)
        js = {
            "run": args.run,
            "cell": cell,
            "organism": organism,
            "celltype": celltype,
            "mids": args.mids,
            "outdir": args.outdir,
            "protocol": args.protocol,
            "barcodes": args.barcodes,
            "umis": args.umis,
            "mismatches": args.mismatches
        }
        js["samples"] = getSamples(myfile)
        writeJson(outfile, js)
