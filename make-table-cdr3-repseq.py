from __future__ import print_function
import sys
import numpy as np

myfile = sys.argv[1]

try:
    fh = open(myfile)
except:
    sys.exit("cannot open file")

reads = list()
percentage = list()
n_samples = 0
for line in fh:
    line = line.strip()
    c = line.split()
    c[2] = c[2].replace("%", "")

    if c[0] == "C":
        n_samples += 1
        if n_samples > 9:  # skip the mouseBCR samples
            print(line)
            reads.append(int(c[1]))
            percentage.append(float(c[2]))

print(myfile)
print("total reads", sum(reads))
print("mean percentage", np.mean(percentage))
print("nr samples", n_samples)
