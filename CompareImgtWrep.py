from __future__ import print_function
import sys


def reformatGeneName(gene):
    '''
    Description: reformat gene name
    In: gene name
    Out: reformatted gene name
    '''

    elem = gene.split(" ")
    genes = list()
    for e in elem:
        if e.startswith("TRB"):
            # remove everything after the asterix and add to list
            name, rest = e.split("*")
            genes.append(name)

    genes.sort()
    genes = list(set(genes))
    newgene = "+".join(genes)

    return(newgene)


def parseFile(f, acc, v, j, cdr3, func, frame):
    '''
    Description: open file en extract info about acc, v, j and cdr3
    In: f (filename) and the indices for the acc,v,j,cdr3 columns
    Out: newfile
    '''
    newfile = f + ".parsed.csv"

    try:
        fh = open(f)
    except:
        sys.exit("cannot open file " + f)

    try:
        fhOut = open(newfile, "w")
        print("acc v j cdr3", file=fhOut)
    except:
        sys.exit("cannot create file " + newfile)

    fh.readline()
    for line in fh:
        line = line.strip()
        c = line.split("\t")
        if func > -1:      # this is an IMGT file
            if c[func] != "productive":
                continue
            if c[frame] != "in-frame":
                continue
            c[acc], ext = c[acc].split("_")
        else:              # this is a WREP file: select cdr3_qual_min >= 30, remove last character from the CDR3
            # if int(c[7]) < 30:
            #     continue
            c[cdr3] = c[cdr3][:-1]
        if " " in c[v]:    # this is an IMGT notation for the V gene
            c[v] = reformatGeneName(c[v])
            c[j] = reformatGeneName(c[j])
        print(c[acc], c[v], c[j], c[cdr3], file=fhOut)

    fh.close()
    fhOut.close()

    return(newfile)


def parseImgt(f):
    '''
    Description: parse file and store the relevant columns
    In: filename
    Out: new file
    '''
    (acc, v, j, cdr3, func, frame) = (1, 3, 9, 20, 2, 21)  # column indices
    newfile = parseFile(f, acc, v, j, cdr3, func, frame)
    return(newfile)


def parseWrep(f):
    '''
    Description: parse file and store the relevant columns
    In: filename
    Out: new file
    '''
    (acc, v, j, cdr3) = (0, 19, 20, 5)   # column indices
    newfile = parseFile(f, acc, v, j, cdr3, -1, -1)
    return(newfile)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit("Usage: " + sys.argv[0] + " 1_Summary.txt rr.all_info.csv")

    imgtFile = sys.argv[1]
    wrepFile = sys.argv[2]

    print(parseImgt(imgtFile))
    print(parseWrep(wrepFile))
