from __future__ import print_function
import sys
import matplotlib.pyplot as plt

# if len(sys.argv) < 2:
#     sys.exit("Usage: plot-xy-clone-overlap.py runXX-clones.csv")


def readSamples(datafile, sampleX, sampleY):

    try:
        fh = open(datafile)
    except:
        sys.exit("cannot load datafile")

    header = fh.readline()
    header = header.strip().split()
    c_sample = header.index("Sample")
    c_v = header.index("V_sub")
    c_j = header.index("J_sub")
    c_cdr3 = header.index("cdr3pep")
    # c_freq = header.index("freq")
    c_perc = header.index("read_perc")

    data = dict()

    # Put all data for sampleX and sampleY in a dictionary
    for line in fh:
        c = line.strip().split()

        if c[c_sample] == sampleX:
            data[(c[c_v], c[c_j], c[c_cdr3])] = data.get((c[c_v], c[c_j], c[c_cdr3]), dict())
            data[(c[c_v], c[c_j], c[c_cdr3])][c[c_sample]] = float(c[c_perc])
        elif c[c_sample] == sampleY:
            data[(c[c_v], c[c_j], c[c_cdr3])] = data.get((c[c_v], c[c_j], c[c_cdr3]), dict())
            data[(c[c_v], c[c_j], c[c_cdr3])][c[c_sample]] = float(c[c_perc])
        else:
            pass

    # Go through dictionary and fill x and y
    x = list()
    y = list()

    for clones in data:
        x.append(data[clones].get(sampleX, 0))
        y.append(data[clones].get(sampleY, 0))

    fh.close()
    return(x, y)


def makeScatter(datafile, sampleX, sampleY):
    plotfile = datafile.split("/")[-1] + "-" + sampleX + "-" + sampleY + ".pdf"
    # plotfile = datafile + "-" + sampleX + "-" + sampleY + ".pdf"

    threshold = 0.5   # threshold HEC
    x, y = readSamples(datafile, sampleX, sampleY)

    fig, ax = plt.subplots()
    ax.scatter(x, y)

    max_xy = max(x + y) + 1
    ax.plot((0, max_xy), (threshold, threshold), '--', color="red")
    ax.plot((threshold, threshold), (0, max_xy), '--', color="red")

    # plt.xlim([-1,max_xy])
    # plt.ylim([-1,max_xy])

    # ax.set_xscale('log')
    # ax.set_yscale('log')

    ax.set_xlabel(sampleX, fontsize=20)
    ax.set_ylabel(sampleY, fontsize=20)
    ax.set_title(datafile)

    ax.grid(True)

    # plt.show()
    try:
        plt.savefig(plotfile)
        print("Wrote", plotfile, "to disk")
    except:
        sys.exit("cannot write plotfile to disk")


makeScatter("/home/barbera/TMP/run06-clones_subs.csv", "S074-166_S173", "S074-188_S185")
