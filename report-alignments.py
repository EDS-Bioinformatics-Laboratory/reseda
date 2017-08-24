from __future__ import print_function
import sys
import os


def countAlignments(sam):
    '''
    Description: counts aligned sequences (unique accessions)
    In: sam (filename)
    Out: count (int)
    '''
    accs = dict()
    try:
        fh = open(sam)
    except:
        sys.exit("cannot open file: " + sam)

    for line in fh:
        line = line.rstrip()
        if line.startswith("@"):  # skip sam header
            continue
        line = line.split("\t")
        if line[1] == "16":       # check if alignment flag is 16 (SEQ being reverse complemented)
            accs[line[0]] = accs.get(line[0], 0) + 1

    count = len(accs)
    return(count)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage:", sys.argv[0], "*.sam")

    fhValign = open("report-ALIGNED-V.txt", "w")
    fhJalign = open("report-ALIGNED-J.txt", "w")
    fhValignAgain = open("report-ALIGNED-AGAIN-V.txt", "w")
    fhJalignAgain = open("report-ALIGNED-AGAIN-J.txt", "w")
    fhOther = open("report-ALIGNED-OTHER.txt", "w")

    for sam in sys.argv[1:]:
        sam_basename = os.path.basename(sam)
        count = countAlignments(sam)
        if "V_" in sam_basename and "-e-clean.sam" in sam_basename:
            print(sam_basename, count, file=fhValign)
        elif "J_" in sam_basename and "-e-clean.sam" in sam_basename:
            print(sam_basename, count, file=fhJalign)
        elif "V_" in sam_basename and "aligned-again" in sam_basename:
            print(sam_basename, count, file=fhValignAgain)
        elif "J_" in sam_basename and "aligned-again" in sam_basename:
            print(sam_basename, count, file=fhJalignAgain)
        else:
            print(sam_basename, count, file=fhOther)

    fhValign.close()
    fhJalign.close()
    fhValignAgain.close()
    fhJalignAgain.close()
