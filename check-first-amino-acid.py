from __future__ import print_function
import sys

if len(sys.argv) < 2:
    sys.exit("Usage: check-first-amino-acid.py *-clones-subs.csv")

for myfile in sys.argv[1:]:
    try:
        fh = open(myfile)
        fhOut = open(myfile+"-cdr3-not-C.txt", "w")

    except:
        sys.exit("cannot open file "+myfile)

    header = fh.readline()
    header = header.strip().split()
    c_cdr3 = header.index("cdr3pep")

    aa = dict()
    total = 0
    for line in fh:
        line = line.strip()
        c = line.split()
        cdr3 = c[c_cdr3][0]
        aa[cdr3] = aa.get(cdr3, 0) + 1
        total += 1
        if cdr3[0] != "C":
            print(line, file=fhOut)

    fh.close()
    fhOut.close()

    for letter in sorted(aa, key=aa.get, reverse=True):
        perc = round(100.0 * aa[letter] / total, 2)
        print(letter, aa[letter], perc)
