from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

# myfile = "/home/barbera/TMP/S074-299-BCRh_S22_L001.assembled-ACTGACTG-IGH_HUMAN-all_info.csv.rr.csv"
# myfile = "/home/barbera/TMP/S074-235_S230_L001.assembled-CGTACGTA-TRB_HUMAN-all_info.csv.rr.csv"
# myfile = "/home/barbera/TMP/S074-002_S73_L001.assembled-ACGTACGT-TRB_HUMAN-all_info.csv.rr.csv"
myfile = "/home/barbera/TMP/S074-228_S223_L001.assembled-ACGTACGT-TRB_HUMAN-all_info.csv.rr.csv"

try:
    fh = open(myfile)
except:
    sys.exit("cannot open file")

top = 10

# Read fractions and put it in list of list
header = fh.readline()
entries = list()        # to store all the data in
for i in range(top):
    entries.append(list())
cdr3s = list()

n = 0
first = list()    # store fractions of V gene that has highest occurrence
for line in fh:
    if n >= 50:
       break

    line = line.strip()
    c = line.split()
    fracs = c[4].split(",")
    fracs = [float(f) for f in fracs]

    # Calculate the cumulative fraction
    for i in range(1,len(fracs)):
        fracs[i] = fracs[i-1] + fracs[i]

    if len(fracs) < top:
        fracs = fracs + (top-len(fracs)) * [1]

    # Store highest occurrence V-gene
    first.append(fracs[0])

    fracs = fracs[0:top]
    fracs.reverse()
    fracs[0] = 1         # first entry should be 1, because that is the "rest" group

    for i in range(top):
        try:
            entries[i].append(fracs[i])
        except:
            entries[i].append(0)

    cdr3s.append(c[0])
    n += 1

fh.close()

# Make the stacked barplot
fig, ax = plt.subplots(figsize=(20, 10)) 
fig.subplots_adjust(bottom=0.3)

N = len(cdr3s)
print("N = ", N)
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

colors = ['gray','black','lightblue','brown','orange','purple','yellow','blue','red','green']
p = list()
for i in range(top):
    data = entries[i]
    p.append(plt.bar(ind, data, width, color=colors[i]))

plt.ylabel('Fraction')
plt.title('V gene assignment')
plt.xticks(ind + width/2., cdr3s, rotation=90)
plt.yticks(np.arange(0, 1.1, 0.1))
labels = ['First', 'Second', 'Third', 'Fourth', 'Fifth', 'Sixth', 'Seventh', 'Eighth', 'Nineth', 'Rest']
labels.reverse()
plt.legend((p[0][0], p[1][0], p[2][0], p[3][0], p[4][0], p[5][0], p[6][0], p[7][0], p[8][0], p[9][0]), labels)
plt.savefig(myfile + ".stackedbars.svg")

# # Now plot histogram of fractions of the most frequent V gene
# fig, ax = plt.subplots(figsize=(20, 10)) 
# hist, bins = np.histogram(first, bins=50)
# width = 0.7 * (bins[1] - bins[0])
# center = (bins[:-1] + bins[1:]) / 2
# plt.bar(center, hist, align='center', width=width)
# plt.xticks(np.arange(0,1.1,0.1))
# plt.savefig(myfile + ".histogram.svg")


# plt.show()
