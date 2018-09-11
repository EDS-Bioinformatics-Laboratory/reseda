from __future__ import print_function
import sys
import json


def AaToImgt(seq):
    '''
    Description: create dictionary for aa position to imgt position
    In: aminoacid sequence
    Out: d[aa_pos] = imgt_pos
    '''
    d = dict()
    aa_pos = 0
    imgt_pos = 1
    for aa in seq:
        if aa != ".":
            d[aa_pos] = (imgt_pos, aa)
            # print(aa, aa_pos, imgt_pos)
            aa_pos += 1
        imgt_pos += 1
    return(d)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: python ImgtNumbering.py ref.table.XX.csv")

    try:
        fh = open(sys.argv[1])
    except:
        sys.exit("cannot open file: " + sys.argv[1])

    header = fh.readline()
    header = header.rstrip()
    header = header.split(",")

    js = dict()
    for line in fh:
        line = line.rstrip()
        c = line.split(",")
        v_gene = c[header.index("V.gene")]
        func = c[header.index("func")]
        seq = c[header.index("seq")]
        js[v_gene] = AaToImgt(seq)

    fh.close()

    # Write new json file
    jsonNew = sys.argv[1].replace(".csv", ".json")
    fhJson = open(jsonNew, "w")
    # print(json.dumps(js, indent=4), file=fhJson)
    print(json.dumps(js), file=fhJson)
    fhJson.close()
