from __future__ import print_function
import sys


def countFirstNucleotides(myfile, nrBases, colNr):
    '''
    Description: open file, get nucleotide sequence from column, count occurrence of first ten nucleotides
    In: filename (str), nrBases (int), column index (int)
    Out: d["ACTGACTG"] = 203894
    '''
    d = dict()
    fh = open(myfile)
    fh.readline()  # skip first line
    for line in fh:
        if "\t" not in line and " " in line:
            raise Warning("input file not tab-delimited")
        c = line.split("\t")
        seq = c[colNr]
        subseq = seq[0:nrBases]
        d[subseq] = d.get(subseq, 0) + 1
    fh.close()
    return(d)


def getHighest(d):
    '''
    Description: return key with highest value
    In: d (dict)
    Out: key, value
    '''
    mykey = sorted(d, key=d.get, reverse=True)[0]
    myvalue = d[mykey]
    return(mykey, myvalue)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " tab-delimited-file-with-sequence (e.g. all_info.csv)")

    for myfile in sys.argv[1:]:
        d = countFirstNucleotides(myfile, 10, 16)
        mid, freq = getHighest(d)
        print(mid, freq, myfile)
