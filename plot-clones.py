from __future__ import print_function
import sys
from bbeeswarm import *
import random

# if len(sys.argv)<2:
#     sys.exit("Usage: plot-clones.py *.subclones.csv")


def remove90PercentSingletons(freqs, percs, colors):
    '''
    Description: remove 90% of the entries with read count 1
    '''
    # First get a list of all indices where the frequency is 1
    indices = list()
    for i in range(len(freqs)):
        if freqs[i] == 1:
            indices.append(i)
    indices = random.sample(indices, int(0.1 * len(indices)))
    indices.sort()

    # Populate new lists
    newfreqs = list()
    newpercs = list()
    newcolors = list()
    for i in range(len(freqs)):
        if i in indices:
            newfreqs.append(freqs[i])
            newpercs.append(percs[i])
            newcolors.append(colors[i])

    return(newfreqs, newpercs, newcolors)


def readData(fileName):
    try:
        fh = open(fileName)
    except:
        sys.exit("cannot open file " + fileName)

    header = fh.readline().strip().split()
    c_freq = header.index("freq")
    c_perc = header.index("read_perc")
    c_vgene = header.index("V_sub")

    freqs = list()
    percs = list()
    colors = list()
    i = 0
    for line in fh:
        # if i > 200:   # temporary for testing
        #     break
        c = line.strip().split()
        freqs.append(int(c[c_freq]))
        percs.append(float(c[c_perc]))
        if c[c_vgene] == "TRBV12-3" or c[c_vgene] == "TRBV12-4":
            colors.append("red")
        elif c[c_vgene] == "TRBV12-3+TRBV12-4":
            colors.append("green")
        else:
            colors.append("black")
        i += 1

    fh.close()

    return(freqs, percs, colors)


imageFile = "beeswarm-before-after-v-reassignment.pdf"
samples = ["S74 before", "S74 after"]
# myfiles = sys.argv[1:]
myfiles = ["/home/barbera/TMP/S074-004_S74_L001.assembled-ACTGACTG-TRB_HUMAN-clones-subs.csv", "/home/barbera/TMP/S074-004_S74_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv.rr.clones_subs.csv"]

allfreqs = list()
allpercs = list()
allcolors = list()

for myfile in myfiles:
    # Read data from the clones file
    (freqs, percs, colors) = readData(myfile)

    # # Make selection, otherwise we need to draw a lot of points
    # (s_freqs, s_percs, s_colors) = remove90PercentSingletons(freqs, percs, colors)  # this one gives an error with KDE calculation
    # (s_freqs, s_percs, s_colors) = (freqs[:3000], percs[:3000], colors[:3000])      # works

    # Add info to the lists that will be provided to bbeeswarm
    allfreqs.append(freqs)
    allpercs.append(percs)
    allcolors.append(colors)

bbeeswarm(samples, allpercs, allcolors, imageFile)
