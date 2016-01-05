"""
Make scatter plot with reads versus umi's
"""
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit("Usage: scatter_reads_umis.py *-clones-subs.csv")

datafiles = sys.argv[1:]

def makeScatter (datafile):
    #plotfile = datafile.split("/")[-1] + ".pdf"
    plotfile = datafile + ".png"

    try:
        data = np.loadtxt(datafile, delimiter='\t', usecols=(3,4), skiprows=1)
    except:
        sys.exit("cannot load datafile with numpy")

    # print(data)

    fig, ax = plt.subplots()
    ax.scatter(data[:,[0]], data[:,[1]])

    ax.set_xlabel("Reads", fontsize=20)
    ax.set_ylabel("UMI's", fontsize=20)
    ax.set_title(datafile)

    ax.grid(True)
    #fig.tight_layout()

    #plt.show()
    try:
        plt.savefig(plotfile)
        print("Wrote", plotfile, "to disk")
    except:
        sys.exit("cannot write plotfile to disk")

############# Main #############

for datafile in datafiles:
    makeScatter(datafile)
