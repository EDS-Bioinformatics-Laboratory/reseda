from __future__ import print_function
import sys
from AmbiguousHlaAlleles import getLookupTable

def parseHlaName (oldname):
    '''
    Description: get relevant part of HLA name
    In: str oldname
    Out: str newname
    '''
    name, alleles = oldname.split("*")
    name = name.split("_")[-1]

    alleles = alleles.split("_")[0]
    alleles = alleles.split(":")
    newname = name + "*" + ":".join(alleles)
    return(newname)

def trimHlaName (oldname):
    '''
    Description: trim HLA name to 6 digits
    In: str oldname
    Out: str newname
    '''
    gene, allele = oldname.split("*")
    alleles = allele.split(":")
    newname = gene + "*" + ":".join(alleles[:3])

    return(newname)

def trimNames (infile, lookupTable):
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

        # Parse relevant bits from the hlaname
        hlaname_new = parseHlaName(hlaname)

        # Lookup the main name of the gene in case it is an ambiguous name
        # If it is not ambiguous use the original name
        hlaname_new = lookupTable.get(hlaname_new, hlaname_new)

        # Trim the HLA name to 6 digits and sum the frequency
        hlaname_new = trimHlaName(hlaname_new)
        newnames[hlaname_new] = newnames.get(hlaname_new, 0) + freq

    for hlaname in sorted(newnames, key=newnames.get, reverse=True):
        print(newnames[hlaname], hlaname, file=fhOut)

    fhIn.close()
    fhOut.close()
    return(outfile)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " *.hla.count.txt")

    lookupTable = getLookupTable("/home/barbera/git/IMGTHLA/xml/hla_ambigs.xml")

    for infile in sys.argv[1:]:
        outfile = trimNames(infile, lookupTable)
        print("Wrote", outfile, "to disk")
