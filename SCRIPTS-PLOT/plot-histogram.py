from __future__ import print_function
import sys
import matplotlib.pyplot as plt
import numpy as np

try:
    fh = open("S074-299-BCRh_S22_L001.assembled-ACTGACTG-IGH_HUMAN-all_info-nr-vgenes.csv")
except:
    sys.exit("cannot open file")

data = list()
header = fh.readline()
for line in fh:
    line = line.strip()
    (cdr3, nr_v_genes, nr_accs) = line.split(",")
    data.append(int(nr_v_genes))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# ax.set_yscale('log')

hist, bins = np.histogram(data, bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.bar(center, hist, align='center', width=width)
plt.show()
