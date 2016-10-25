from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

def bubblePlot (f,d):
    '''
    Description: Make bubble plot
    In: f (filename), d[(a,b)] = size, where a and b are names (e.g. gene names)
    Out: plot
    '''

    # Get the unique names for x and y axis
    names = list()
    for a,b in d:
        names.append(a)
        names.append(b)
    names = list(set(names))
    names.sort()

    # Fill lists with x, y and size values
    x = list()
    y = list()
    s = list()
    for a, b in d:
        x.append(names.index(a))
        y.append(names.index(b))
        s.append(d[(a, b)])

    # Define colors
    cmap = plt.cm.get_cmap('cool')
    ncolors = max(s) + 1
    colors = [cmap(i*(256/ncolors)) for i in range(ncolors)]
    usecolors = [ colors[size] for size in s ]

    # Make figure
    fig, ax = plt.subplots(figsize=(20, 20)) 
    # fig.subplots_adjust(bottom=0.3)
    plt.scatter(x,y,s=s,c=usecolors)
    plt.xticks(range(len(names)), names, rotation=90)
    plt.yticks(range(len(names)), names)
    image = f + ".bubble.svg"
    plt.savefig(image)
    return(image)

if __name__ == '__main__':
    myfile = "/home/barbera/Data/tbcell/RepSeq2016/IMGT-pairwise/TRBV_human.distances.txt"
    try:
        fh = open(myfile)
    except:
        sys.exit("cannot open file: " + myfile)

    d = dict()

    header = fh.readline()
    for line in fh:
        line = line.strip()
        (a, b, nope, size) = line.split()

        # a_sub, alleleA = a.split("*")
        # b_sub, alleleB = b.split("*")
        # if a_sub == b_sub:  # Skip same gene groups
        #     continue

        d[(a,b)] = int(size)
        d[(b,a)] = int(size)

    fh.close()

    print("Wrote", bubblePlot(myfile,d), "to disk")
