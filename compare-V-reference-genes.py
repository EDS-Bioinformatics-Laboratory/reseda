from __future__ import print_function
import sys
from PairwiseAlignment import pairAlign

def alignAllVsAll (fileOut,sequences):
    fhOut = open(fileOut, "w")
    seqs = sequences.keys()   # [0:10]
    for i in range(0,len(seqs)-1):
        idA = seqs[i]
        seqA = sequences[idA]

        for j in range(i+1,len(seqs)):
            idB = seqs[j]
            seqB = sequences[idB]
            # print(idA, idB, seqA, seqB)
            print(idA, idB)
            aln = pairAlign(seqA, seqB)
            print(idA, idB, aln.distance, aln.first, aln.second, file=fhOut)
    fhOut.close()

def compareGenes (reftable):
    try:
        fh = open(reftable)
    except:
        sys.exit("cannot open file " + reftable)

    vgenes = dict()

    # Check the column names, get index of the columns "V.gene"/"name" and "seq"
    header = fh.readline()
    header = header.replace('"','')
    header = header.rstrip().split(",")
    if "V.gene" in header:
        c_name = header.index("V.gene")
    elif "name" in header:
        c_name = header.index("name")
    else:
        sys.exit("ERR: unknown which column contains the gene name")
    if "seq" in header:
        c_seq = header.index("seq")
    else:
        sys.exit("ERR: unknown which column contains the sequence")

    # Read the sequences
    for line in fh:
        line = line.replace('"','')
        line = line.rstrip()
        c = line.split(",")
        name = c[c_name]
        seq = c[c_seq]
        vgenes[name] = seq

    fh.close()

    fileOut = reftable + ".pairalign.txt"
    alignAllVsAll(fileOut, vgenes)

########## MAIN ############

reftables = ["ref.table.BCRk.csv","ref.table.BCRl.csv","ref.table.heavy.csv","ref.table.Ig.csv","ref.table.mouse.BCRk.csv","ref.table.mouse.BCRl.csv","ref.table.mouse.heavy.csv","ref.table.TCRa.csv","ref.table.TCRb.csv"]
for reftable in reftables:
    compareGenes(reftable)
