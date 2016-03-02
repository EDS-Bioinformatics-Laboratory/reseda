from __future__ import print_function
import sys

if len(sys.argv) < 3:
    sys.exit("Usage: select-gene-from-sam.py genename samfile(s)")

def selectGene (gene, sam):
    fh = open(sam, "r")
    outfile = sam+"-"+gene+".sam"
    fhOut = open(outfile, "w")
    for line in fh:
        line = line.rstrip()
        if line.startswith("@"):
            continue
        if gene in line:
            print(line, file=fhOut)
    fh.close()
    fhOut.close()
    print("wrote", outfile, "to disk")

########### MAIN ############

gene = sys.argv[1]

for myfile in sys.argv[2:]:
    selectGene(gene, myfile)
