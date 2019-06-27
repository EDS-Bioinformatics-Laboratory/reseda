from __future__ import print_function
import sys
import matplotlib.pyplot as plt
import numpy as np


def readData(datafile, colname, DELIMITER):
    '''
    Description: read samplename and cdr3 from file and store in dictionary
    In: filename, column with cdr3, field delimiter
    Out: data[(sample,cdr3)] = count (count should be one)
    '''
    try:
        fh = open(datafile)
    except:
        sys.exit("cannot open file: " + datafile)

    data = dict()
    header = fh.readline().strip()
    header = header.split(DELIMITER)
    c_get = header.index(colname)
    c_sample = header.index("Sample")

    for line in fh:
        line = line.strip()
        c = line.split(DELIMITER)
        if "BCRh" in c[c_sample]:   # only take the BCRh or TCRb samples (BCRh in or not in c[c_sample])
            data[(c[c_sample], c[c_get])] = data.get((c[c_sample], c[c_get]), 0) + 1
            # if data[(c[c_sample],c[c_get])] > 1: # check if we have seen this combination before. If so, print a WARNING
            #     print("WARNING:", datafile, c[c_sample], c[c_get], "is not unique!")

    fh.close()
    return(data)


def removeReassignedClones(d1, d2):
    '''
    Description: Remove the reassigned clones (d2) from d1
    In: d1[sample,cdr3], d2[sample,cdr3]
    Out: d1[sample,cdr3]
    '''
    count_discarded = 0
    d2_new = dict()  # keep only those cdr3's that were not filtered out due to low quality CDR3
    for (sample, cdr3) in d2:
        if (sample, cdr3) in d1:
            del d1[(sample, cdr3)]  # delete sample and cdr3 from the clones_subs dictionary
            d2_new[(sample, cdr3)] = d2[(sample, cdr3)]  # keep in reassigned dictionary
        else:
            print("WARNING:", sample, cdr3, "not in clones_sub")
            count_discarded += 1
    print("WARNING:", count_discarded, "sample,cdr3 combinations removed from reassigned dictionary because these were not present in clones_sub")

    return(d1, d2_new)


def getLength(d):
    '''
    Description: check length of cdr3
    In: d[(sample,cdr3)]
    Out: lengths (list with integers)
    '''

    data = list()

    for (sample, cdr3) in d:
        length = len(cdr3)
        data.append(length)

    return(data)


if __name__ == '__main__':

    if len(sys.argv) < 3:
        # e.g. run06-clones_subs.csv run06.rr.csv
        sys.exit("Usage: " + sys.argv[0] + " SAMPLE.rr.clones_subs.csv SAMPLE.rr.csv")

    all_clones = sys.argv[1]
    reassigned_clones = sys.argv[2]

    # Read files and store: sample, cdr3 (LET OP: er wordt gefilterd obv de sample naam in deze functies!)
    dict1 = readData(all_clones, "cdr3pep", "\t")
    dict2 = readData(reassigned_clones, "cdr3", "\t")

    # Remove all entries from dict1 that are in dict2 (remove the reassigned clones)
    (dict1, dict2) = removeReassignedClones(dict1, dict2)

    # Get CDR3 lengths
    d1 = getLength(dict1)
    d2 = getLength(dict2)

    # Write CDR3 lengths to a file
    fhAll = open("cdr3-length-all-but-reassigned.txt", "w")
    fhReassigned = open("cdr3-length-reassigned.txt", "w")
    for d in d1:
        print(d, file=fhAll)
    for d in d2:
        print(d, file=fhReassigned)
    fhAll.close()
    fhReassigned.close()

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_yscale('log')

    binBoundaries = np.linspace(0, 100, 101)

    plt.hist(d1, bins=binBoundaries, color="black", label="All clones")
    plt.hist(d2, bins=binBoundaries, color="yellow", label="Clones with re-assigned V", hatch="\\")

    plt.title("CDR3 length distribution")
    plt.xlabel("CDR3 peptide length")
    plt.ylabel("Frequency")
    plt.legend()
    plt.show()
