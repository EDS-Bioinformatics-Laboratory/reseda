from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

pear = "report-PEAR.txt"
mids = "report-MIDs.txt"
cdr3 = "report-CDR3.txt"
prod = "report-PRODUCTIVE.txt"
reassign = "report-AFTER-V-REASSIGNMENT.txt"

def parsePear (pear):
    try:
        fh = open(pear)
    except:
        sys.exit("cannot open file " + pear)

    totalreads = dict()
    summary = dict()
    for line in fh:
        line = line.strip()
        c = line.split()
        path = c[0].split("/")[-1]
        sample,rest = path.split("_L001")
        total = int(c[-2].replace(",",""))
        assembled = int(c[-4].replace(",",""))
        percentage = 100.0 * assembled / total

        # print(sample, assembled, total, round(percentage, 2))
        totalreads[sample] = total
        summary[sample] = (assembled, round(percentage, 2))
    return(totalreads, summary)

def parseMids (totalreads, mids):
    try:
        fh = open(mids)
    except:
        sys.exit("cannot open file " + mids)

    summary = dict()
    midcount = dict()
    for line in fh:
        line = line.strip()
        if line.startswith(":"):   # skip these lines
            continue
        elif line.startswith("/mnt"): # parse sample name
            line = line.split("/")[-1]
            sample, rest = line.split("_L001")
            midcount[sample] = dict()
        else:                      # parse MIDs and frequency
            mid, freq = line.split()
            if mid != "nomatch":
                midcount[sample][mid] = int(freq)

    for sample in midcount:
        for mid in sorted(midcount[sample], key=midcount[sample].get, reverse=True):
            percentage = round(100.0 * midcount[sample][mid] / totalreads[sample],2)
            summary[sample] = (mid, midcount[sample][mid], percentage)
            break

    return(summary)

def parseCdr3 (totalreads, summary_mids, cdr3):
    try:
        fh = open(cdr3)
    except:
        sys.exit("cannot open file " + cdr3)

    summary = dict()
    for line in fh:
        line = line.strip()
        c = line.split()
        path = c[0].split("/")[-1]
        sample, rest = path.split("_L001.assembled-")
        mid = rest.split(".")[0]
        freq = int(c[-2])
        if mid in summary_mids[sample]:
            percentage = round(100.0 * freq / totalreads[sample], 2)
            summary[sample] = (freq, percentage)

    return(summary)

def parseProductive (totalreads, summary_mids, prod):
    try:
        fh = open(prod)
    except:
        sys.exit("cannot open file " + prod)

    summary = dict()
    for line in fh:
        line = line.strip()
        c = line.split()
        path = c[0].split("/")[-1]
        sample, rest = path.split("_L001.assembled-")
        mid = rest.split("-")[0]
        freq = int(c[-1])
        if mid in summary_mids[sample]:
            percentage = round(100.0 * freq / totalreads[sample], 2)
            summary[sample] = (freq, percentage)

    return(summary)

def parseReassign (totalreads, summary_mids, reassign):
    try:
        fh = open(reassign)
    except:
        sys.exit("cannot open file " + reassign)

    summary = dict()
    for line in fh:
        line = line.strip()
        c = line.split()
        path = c[0].split("/")[-1]
        sample, rest = path.split("_L001.assembled-")
        mid = rest.split("-")[0]
        freq = int(c[-1])
        if mid in summary_mids[sample]:
            percentage = round(100.0 * freq / totalreads[sample], 2)
            summary[sample] = (freq, percentage)

    return(summary)

def makeBarChart (plotfile,title,y_label, threshold,x,y,z,a,b):
    x_pos = np.arange(len(x))
    
    fig, ax = plt.subplots(figsize=(60, 10)) 
    fig.subplots_adjust(bottom=0.3)
    p = list()
    p.append(ax.bar(x_pos, y, align='center', color="red"))
    p.append(ax.bar(x_pos, z, align='center', color="blue"))
    p.append(ax.bar(x_pos, a, align='center', color="yellow"))
    p.append(ax.bar(x_pos, b, align='center', color="green"))
    ax.plot((0,max(x_pos)),(threshold,threshold), '--', color="black")

    plt.xticks(x_pos, x, rotation=90)
    ax.set_xlabel('Samples')
    ax.set_ylabel(y_label)
    ax.set_title(title)
    labels = ['With correct MID', 'CDR3 identified', 'VJ assigned', 'After V re-assignment']
    plt.legend((p[0][0], p[1][0], p[2][0], p[3][0]), labels)

    
    try:
        fig.savefig(plotfile)
        print("Wrote", plotfile, "to disk")
    except:
        sys.exit("cannot write plotfile to disk")

############# MAIN ################

# Parse log files
(totalreads, summary_pear) = parsePear(pear)
summary_mids = parseMids(totalreads, mids)
summary_cdr3 = parseCdr3(totalreads, summary_mids, cdr3)
summary_prod = parseProductive(totalreads, summary_mids, prod)
summary_reassign = parseReassign(totalreads, summary_mids, reassign)

# Create one big table with all summary statistics
try:
    fhOut = open("report-all.csv", "w")
except:
    sys.exit("cannot write file report-all.csv")

print("Sample TotalReads MID MidFreq MidPerc Cdr3Freq Cdr3Perc VJFreq VJPerc ReassignedFreq ReassignedPerc", file=fhOut)
samples = list()
midpercs = list()
cdr3percs = list()
prodpercs = list()
reassigns = list()
midfreqs = list()
cdr3freqs = list()
prodfreqs = list()
reassignfreqs = list()
for sample in sorted(totalreads):
    # Print all numbers to file
    total = totalreads[sample]
    (mid, midfreq, midperc) = summary_mids[sample]
    (cdr3freq, cdr3perc) = summary_cdr3.get(sample,(0,0))
    (prodfreq, prodperc) = summary_prod.get(sample,(0,0))
    (reassignfreq, reassignperc) = summary_reassign.get(sample,(0,0))
    print(sample, total, mid, midfreq, midperc, cdr3freq, cdr3perc, prodfreq, prodperc, reassignfreq, reassignperc, file=fhOut)

    # Store percentages in lists
    samples.append(sample)
    midpercs.append(midperc)
    cdr3percs.append(cdr3perc)
    prodpercs.append(prodperc)
    reassigns.append(reassignperc)

    # Store frequencies in lists
    midfreqs.append(midfreq)
    cdr3freqs.append(cdr3freq)
    prodfreqs.append(prodfreq)
    reassignfreqs.append(reassignfreq)

makeBarChart("report-all-percentages.pdf", "Run06", "Reads (percentage)", 70, samples,midpercs,cdr3percs,prodpercs,reassigns)
makeBarChart("report-all-frequencies.pdf", "Run06", "Reads (frequency)", 50000, samples,midfreqs,cdr3freqs,prodfreqs,reassignfreqs)
