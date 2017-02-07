from __future__ import print_function
import sys
import random

if len(sys.argv) < 3:
    sys.exit("Usage: split-data-at-random.py 2 clones.csv")

for i in range(10):
    try:
        nr_samples = int(sys.argv[1])
    except:
        sys.exit("Usage: split-data-at-random.py 2 clones.csv")

    try:
        fh = open(sys.argv[2])
    except:
        sys.exit("cannot open file")

    try:
        fhOut = open("clones" + str(i) + ".csv", "w")
    except:
        sys.exit("cannot open file")

    header = fh.readline()
    for line in fh:
        line = line.replace('"', '')
        line = line.rstrip()
        c = line.split(",")
        # v = c[0]
        # j = c[1]
        # cdr3 = c[2]
        # vjcdr3 = "|".join(c[0:3])
        vjcdr3 = c[8]
        freq = int(c[9])

        # new frequencies
        new_freq = [0] * nr_samples

        # divide frequency over nr_samples
        for knikker in range(freq):
            dice = random.randint(0, nr_samples - 1)
            new_freq[dice] += 1

        print(",".join(c + [str(c) for c in new_freq]), file=fhOut)
    fh.close()
    fhOut.close()
