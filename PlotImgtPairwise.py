from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

def reformatData(d):
    '''
    Description: Reformat data for use in a figure
    In: d[(a,b)]
    Out: x, y, s
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

    # Set diagonal to a high value
    for i in range(len(names)):
        x.append(i)
        y.append(i)
        s.append(200)

    # Define colors
    cmap = plt.cm.get_cmap('cool')
    ncolors = max(s) + 1
    colors = [cmap(i*(256/ncolors)) for i in range(ncolors)]
    usecolors = [ colors[size] for size in s ]

    return(x, y, s, usecolors, names)

def asArray (x, y, s, names):
    '''
    Description: Data in numpy array format
    In: x, y, s
    Out: X, Y, C
    '''

    C = np.zeros((len(names), len(names)))
    for i in range(len(x)):
        C[x[i]][y[i]] = s[i]
    C = np.array(C)

    return(C)

def bubblePlot (f, x, y, s, usecolors, names):
    '''
    Description: Make bubble plot
    In: f (filename), x, y, s, names
    Out: plot
    '''

    # Make figure
    fig, ax = plt.subplots(figsize=(20, 20)) 
    plt.scatter(x,y,s=s,c=usecolors)
    plt.xticks(range(len(names)), names, rotation=90)
    plt.yticks(range(len(names)), names)
    image = f + ".bubble.svg"
    plt.savefig(image)
    return(image)

def heatMap (f, C, names):
    '''
    Description: make heatmap
    In: f (filename), d[(a,b)] = value
    Out: image
    Code adapted from: https://nbviewer.jupyter.org/gist/joelotz/5427209
    '''
    # Make figure
    fig, ax = plt.subplots(figsize=(20, 20)) 
    # plt.pcolor(X, Y, C, cmap='cool')
    plt.pcolor(C, cmap='hot')

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(C.shape[0])+0.5, minor=False)
    ax.set_xticks(np.arange(C.shape[1])+0.5, minor=False)

    # want a more natural, table-like display
    # ax.invert_yaxis()
    # ax.xaxis.tick_top()

    # Set labels
    ax.set_xticklabels(names, minor=False) 
    ax.set_yticklabels(names, minor=False)

    # rotate the 
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks(): 
        t.tick1On = False 
        t.tick2On = False 
    for t in ax.yaxis.get_major_ticks(): 
        t.tick1On = False 
        t.tick2On = False  

    # Add legend
    plt.colorbar()

    # Write to file
    image = f + ".heatmap.svg"
    plt.savefig(image)
    return(image)

if __name__ == '__main__':
    # mydir = "/home/barbera/Data/tbcell/RepSeq2016/IMGT-pairwise/"
    mydir = "/home/narya/Werk/RepSeq2016/"
    myfile = mydir + "TRBV_human.distances.txt"

    try:
        fh = open(myfile)
    except:
        sys.exit("cannot open file: " + myfile)

    d = dict()

    header = fh.readline()
    for line in fh:
        line = line.strip()
        (a, b, size, nope) = line.split()

        a_sub, alleleA = a.split("*")
        b_sub, alleleB = b.split("*")
        if a_sub == b_sub:  # Skip same gene groups
            size = 200

        d[(a,b)] = int(size)
        d[(b,a)] = int(size)

    fh.close()

    (x, y, s, usecolors, names) = reformatData(d)
    C = asArray(x,y,s,names)

    # print("Wrote", bubblePlot(myfile, x, y, s, usecolors, names), "to disk")
    print("Wrote", heatMap(myfile, C, names), "to disk")
