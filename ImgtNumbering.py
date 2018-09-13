from __future__ import print_function
import sys
import json
import argparse


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
    parser = argparse.ArgumentParser(description='Create a json file with translation table for AA position to IMGT numbering from ref.table')
    parser.add_argument('-c', '--chain', default='IGH', type=str, help='IGH, IGK or IGL (default: %(default)s)')
    parser.add_argument('-i', '--input_file', default='ref.table.XX.csv', type=str, help='reference table with IMGT AA sequences with gaps (default: %(default)s)')
    args = parser.parse_args()

    if args.input_file == 'ref.table.XX.csv' or args.chain not in ['IGH', 'IGK', 'IGL']:
        parser.print_help()
        exit()

    print("Chain:", args.chain)
    print("Input:", args.input_file)

    try:
        fh = open(args.input_file)
    except:
        sys.exit("cannot open file: " + args.input_file)

    header = fh.readline()
    header = header.rstrip()
    header = header.replace('"', '')
    header = header.split(",")
    # print("header:", header)

    js = dict()
    for line in fh:
        line = line.rstrip()
        line = line.replace('"', '')
        c = line.split(",")

        # Check if the gene or V.gene column is present: store the v_gene name
        if "V.gene" in header:
            v_gene = c[header.index("V.gene")]
        elif "name" in header:
            v_gene = c[header.index("name")]
        else:
            sys.exit("gene column not found in: " + args.input_file)

        # Convert the V gene name for IGK and IGL
        if args.chain in ['IGK', 'IGL']:
            v_gene = args.chain + "V" + v_gene.replace('.', '-')

        seq = c[header.index("seq")]
        js[v_gene] = AaToImgt(seq)

    fh.close()

    # Write new json file
    jsonNew = args.input_file.replace(".csv", ".json")
    fhJson = open(jsonNew, "w")
    # print(json.dumps(js, indent=4), file=fhJson)
    print(json.dumps(js), file=fhJson)
    fhJson.close()
    print("Wrote", jsonNew, "to disk")
