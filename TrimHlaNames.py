from __future__ import print_function
import sys

def trimNames (infile):
    '''
    Description: trim the HLA names to 6 digits
    In: infile (filename)
    Out: outfile (filename)
    '''
    outfile = infile + ".trimmed.txt"

    try:
        fhIn = open(infile)
        fhOut = open(outfile, "w")
    except:
        sys.exit("can't read or write file: " + infile + ", " + outfile)

    newnames = dict()
    for line in fhIn:
        line = line.strip()  # remove leading whitespace and newline
        (freq, hlaname) = line.split()
        freq = int(freq)

        name, alleles = hlaname.split("*")
        name = name.split("_")[-1]

        alleles = alleles.split("_")[0]
        alleles = alleles.split(":")
        hlaname_new = name + "*" + ":".join(alleles[:3])

        newnames[hlaname_new] = newnames.get(hlaname_new, 0) + freq

    for hlaname in sorted(newnames, key=newnames.get, reverse=True):
        print(newnames[hlaname], hlaname, file=fhOut)

    fhIn.close()
    fhOut.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " *.hla.count.txt")

    for infile in sys.argv[1:]:
        outfile = trimNames(infile)