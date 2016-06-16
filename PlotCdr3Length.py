from __future__ import print_function
import sys
import matplotlib.pyplot as plt
import numpy as np

def getLength (datafile, colname, DELIMITER):
    '''
    Description: check length of content of column specified in "colname" in file "datafile". Delimiter usually "\t", " ", or ","
    In: filename (str), column name (str), DELIMITER (str)
    Out: lengths (list with integers)
    '''
    try:
        fh = open(datafile)
    except:
        sys.exit("cannot open file: " + datafile)

    data = list()
    header = fh.readline().strip()
    header = header.split(DELIMITER)
    c_get = header.index(colname)

    for line in fh:
        line = line.strip()
        c = line.split(DELIMITER)
        length = len(c[c_get])
        data.append(length)

    return(data)

if __name__ == '__main__':

    if len(sys.argv) < 3:
        sys.exit("Usage: " + sys.argv[0] + " SAMPLE.rr.clones_subs.csv SAMPLE.rr.csv")

    all_clones = sys.argv[1]
    reassigned_clones = sys.argv[2]

    d1 = getLength(all_clones, "cdr3pep", "\t")
    d2 = getLength(reassigned_clones, "cdr3", "\t")

    # plt.hist(d1)
    # plt.hist(d2)
    # plt.title("CDR3 length distribution")
    # plt.xlabel("CDR3 length")
    # plt.ylabel("Frequency")
    # plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_yscale('log')

    hist, bins = np.histogram(d1, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width, color="black", label="All clones")

    hist, bins = np.histogram(d2, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width, color="yellow", label="Clones with re-assigned V")

    plt.title("CDR3 length distribution")
    plt.xlabel("CDR3 peptide length")
    plt.ylabel("Frequency")
    plt.legend()
    plt.show()
