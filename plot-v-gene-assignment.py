from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

# myfile = "/home/barbera/TMP/S074-299-BCRh_S22_L001.assembled-ACTGACTG-IGH_HUMAN-all_info.csv.rr.csv"
# myfile = "/home/barbera/TMP/S074-235_S230_L001.assembled-CGTACGTA-TRB_HUMAN-all_info.csv.rr.csv"
# myfile = "/home/barbera/TMP/S074-002_S73_L001.assembled-ACGTACGT-TRB_HUMAN-all_info.csv.rr.csv"
# myfile = "/home/barbera/TMP/S074-228_S223_L001.assembled-ACGTACGT-TRB_HUMAN-all_info.csv.rr.csv"

# myfile = "/home/barbera/TMP/S074-004_S74_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv.rr.csv"  # voor paper
myfile = "/home/barbera/TMP/S074-241-BCRh_S17_L001.assembled-CAGTCAGT-IGH_HUMAN-all_info.csv.rr.csv"  # voor paper

# myfile = sys.argv[1]

try:
    fh = open(myfile)
except:
    sys.exit("cannot open file")

top = 10  # show count for top X V genes, anything above that will be counted as "rest"


def autolabel(rects, labels):
    # attach some text labels
    i = 0
    for rect in rects:
        height = rect.get_height()
        # ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, '%d' % int(height), ha='center', va='bottom')
        ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, '%d' % labels[i], ha='center', va='bottom', rotation=90)
        i += 1


# Read fractions and put it in list of list
header = fh.readline()
entries = list()        # to store all the data in
for i in range(top):
    entries.append(list())
cdr3s = list()
freqs = list()
vgenes_cooccurrence = dict()
all_vgenes = list()

n = 0
first = list()    # store fractions of V gene that has highest occurrence
for line in fh:

    line = line.strip()
    c = line.split()
    fracs = c[4].split(",")
    fracs = [float(f) for f in fracs]

    # Which vgenes occur together
    vgenes = c[2].split(",")
    all_vgenes += vgenes
    for a in range(len(vgenes)):
        for b in range(len(vgenes)):
            if vgenes[a] == vgenes[b]:
                continue
            vgenes_cooccurrence[(vgenes[a], vgenes[b])] = vgenes_cooccurrence.get((vgenes[a], vgenes[b]), 1) + 1

    if n < 50:
        # Calculate the cumulative fraction
        for i in range(1, len(fracs)):
            fracs[i] = fracs[i - 1] + fracs[i]

        if len(fracs) < top:
            fracs = fracs + (top - len(fracs)) * [1]

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
        freqs_tmp = c[3].split(",")
        freqs_tmp = [int(freq) for freq in freqs_tmp]
        freqs.append(sum(freqs_tmp))

    n += 1

fh.close()

# Make the stacked barplot
fig, ax = plt.subplots(figsize=(20, 10))
fig.subplots_adjust(bottom=0.3)

N = len(cdr3s)
print("N = ", N)
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence

colors = ['gray', 'black', 'lightblue', 'brown', 'orange', 'purple', 'yellow', 'blue', 'red', 'green']
p = list()
for i in range(top):
    data = entries[i]
    p.append(plt.bar(ind, data, width, color=colors[i]))

plt.ylabel('Fraction')
plt.title('V gene assignment')
plt.xticks(ind + width / 2., cdr3s, rotation=90)
plt.yticks(np.arange(0, 1.1, 0.1))
labels = ['First', 'Second', 'Third', 'Fourth', 'Fifth', 'Sixth', 'Seventh', 'Eighth', 'Nineth', 'Rest']
labels.reverse()
plt.legend((p[0][0], p[1][0], p[2][0], p[3][0], p[4][0], p[5][0], p[6][0], p[7][0], p[8][0], p[9][0]), labels)
autolabel(p[0], freqs)

plt.savefig(myfile + ".stackedbars.svg")

# Make scatterplot of co-occurring V genes, size is nr of times they occur together
x = list()
y = list()
s = list()
vs = list(set(all_vgenes))
vs.sort()
for v_a, v_b in vgenes_cooccurrence:
    x.append(vs.index(v_a))
    y.append(vs.index(v_b))
    s.append(vgenes_cooccurrence[(v_a, v_b)])

fig, ax = plt.subplots(figsize=(20, 20))
# fig.subplots_adjust(bottom=0.3)
plt.scatter(x, y, s=s)
plt.xticks(range(len(vs)), vs, rotation=90)
plt.yticks(range(len(vs)), vs)
plt.savefig(myfile + ".co-occurrence.svg")
