from __future__ import print_function
import sys
import csv

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage:", sys.argv[0], "pt.table.csv aa.reads.csv")
        exit()

    # Read the pt.table
    fhOut = dict()
    with open(sys.argv[1]) as csvfile:
        fh = csv.reader(csvfile)
        row = 0
        for c in fh:
            if row == 0:
                header = c
            else:
                pt = c[header.index("pt")]
                group = c[header.index("group")]
                group2 = group.replace(" ", "")
                fhOut[(pt, group)] = open(pt + "-" + group2 + "_S" + str(row) + "_L001.assembled.txt", "w")
            row += 1

    # Open the AA.reads file and write entries to the correct output file
    with open(sys.argv[2]) as csvfile:
        fh = csv.reader(csvfile)
        row = 0
        for c in fh:
            if row == 0:
                header = c
            else:
                pt = c[header.index("pt")]
                group = c[header.index("group")]
                acc = c[header.index("accession")]
                seq = c[header.index("seq")]
                print("\t".join([acc, seq]), file=fhOut[(pt, group)])
            row += 1

    for myindex, myfile in fhOut.items():
        myfile.close()
