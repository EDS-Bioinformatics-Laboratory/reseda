from __future__ import print_function
import sys
import json
import numpy as np
import matplotlib.pyplot as plt

pear = "report-PEAR.txt"
mids = "report-MIDs.txt"
cdr3 = "report-CDR3.txt"
prod = "report-PRODUCTIVE.txt"
reassign = "report-AFTER-V-REASSIGNMENT.txt"


def parsePear(pear):
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
        sample, rest = path.split("_L001")
        total = int(c[-2].replace(",", ""))
        assembled = int(c[-4].replace(",", ""))
        percentage = 100.0 * assembled / total

        # print(sample, assembled, total, round(percentage, 2))
        totalreads[sample] = total
        summary[sample] = (assembled, round(percentage, 2))
    return(totalreads, summary)


def parseMids(totalreads, mids):
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
        elif line.startswith("/mnt"):  # parse sample name
            line = line.split("/")[-1]
            sample, rest = line.split("_L001")
            midcount[sample] = dict()
        else:                      # parse MIDs and frequency
            mid, freq = line.split()
            if mid != "nomatch":
                midcount[sample][mid] = int(freq)

    for sample in midcount:
        for mid in sorted(midcount[sample], key=midcount[sample].get, reverse=True):
            percentage = round(100.0 * midcount[sample][mid] / totalreads.get(sample, -1), 2)
            summary[sample] = (mid, midcount[sample][mid], percentage)
            break

    return(summary)


def parseCdr3(totalreads, summary_mids, cdr3):
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
            percentage = round(100.0 * freq / totalreads.get(sample, -1), 2)
            summary[sample] = (freq, percentage)

    return(summary)


def parseProductive(totalreads, summary_mids, prod):
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
            percentage = round(100.0 * freq / totalreads.get(sample, -1), 2)
            summary[sample] = (freq, percentage)

    return(summary)


def parseReassign(totalreads, summary_mids, reassign):
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
            percentage = round(100.0 * freq / totalreads.get(sample, -1), 2)
            summary[sample] = (freq, percentage)

    return(summary)


def makeBarChart(plotfile, title, y_label, threshold, x, v, w, y, z, a, b):
    x_pos = np.arange(len(x))

    fig, ax = plt.subplots(figsize=(60, 10))
    fig.subplots_adjust(bottom=0.3)
    cmap = plt.cm.get_cmap('YlGn')
    ncolors = 5
    colors = [cmap(i * (256 / ncolors)) for i in range(ncolors)]

    p = list()
    p.append(ax.bar(x_pos, v, align='center', color=colors[0], label='Total'))
    p.append(ax.bar(x_pos, w, align='center', color=colors[1], label='Assembled (PEAR)'))
    p.append(ax.bar(x_pos, y, align='center', color=colors[2], label='With correct MID'))
    p.append(ax.bar(x_pos, z, align='center', color=colors[3], label='CDR3 identified'))
    # p.append(ax.bar(x_pos, a, align='center', color=colors[4], label='VJ assigned'))
    p.append(ax.bar(x_pos, b, align='center', color=colors[4], label='VJ assigned'))
    ax.plot((0, max(x_pos)), (threshold, threshold), '--', color="black")

    plt.xticks(x_pos, x, rotation=90)
    ax.set_xlabel('Samples')
    ax.set_ylabel(y_label)
    ax.set_title(title)
    plt.legend()

    try:
        fig.savefig(plotfile)
        print("Wrote", plotfile, "to disk")
    except:
        sys.exit("cannot write plotfile to disk")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage:", sys.argv[0], "2017etcetc-new.json (json file with run and sample information)")
        exit()

    jsonFile = sys.argv[1]  # "20170224_Rheuma_MiSeqRUN013.json"

    # Read json with parsed sample sheet info (made with MetaData.py)
    try:
        fh = open(jsonFile)
    except:
        sys.exit("cannot open file: " + jsonFile)
    text = fh.read()
    js = json.loads(text)
    fh.close()

    # Get all expected sample names
    expected_samples = list()
    for sample in js["Samples"]:
        expected_samples.append(sample["Sample_Name"] + "_" + sample.get("Sample_Nr", ""))
    expected_samples.sort()

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

    print("Sample TotalReads AssembledFreq AssembledPerc MID MidFreq MidPerc Cdr3Freq Cdr3Perc VJFreq VJPerc ReassignedFreq ReassignedPerc", file=fhOut)
    samples = list()
    totals = list()
    assembledfreqs = list()
    assembledpercs = list()
    midpercs = list()
    cdr3percs = list()
    prodpercs = list()
    reassigns = list()
    midfreqs = list()
    cdr3freqs = list()
    prodfreqs = list()
    reassignfreqs = list()
    for sample in expected_samples:
        # Print all numbers to file
        total = totalreads.get(sample, -1)
        (assembledfreq, assembledperc) = summary_pear.get(sample, (0, 0))
        (mid, midfreq, midperc) = summary_mids.get(sample, ("ABC", 0, 0))
        (cdr3freq, cdr3perc) = summary_cdr3.get(sample, (0, 0))
        (prodfreq, prodperc) = summary_prod.get(sample, (0, 0))
        (reassignfreq, reassignperc) = summary_reassign.get(sample, (0, 0))
        print(sample, total, assembledfreq, assembledperc, mid, midfreq, midperc, cdr3freq, cdr3perc, prodfreq, prodperc, reassignfreq, reassignperc, file=fhOut)

        # Store percentages in lists
        samples.append(sample)
        assembledpercs.append(assembledperc)
        midpercs.append(midperc)
        cdr3percs.append(cdr3perc)
        prodpercs.append(prodperc)
        reassigns.append(reassignperc)

        # Store frequencies in lists
        totals.append(total)
        assembledfreqs.append(assembledfreq)
        midfreqs.append(midfreq)
        cdr3freqs.append(cdr3freq)
        prodfreqs.append(prodfreq)
        reassignfreqs.append(reassignfreq)

    makeBarChart("report-all-percentages.pdf", "Run overview", "Reads (percentage)", 70, samples, len(samples) * [100], assembledpercs, midpercs, cdr3percs, prodpercs, reassigns)
    makeBarChart("report-all-frequencies.pdf", "Run overview", "Reads (frequency)", 0, samples, totals, assembledfreqs, midfreqs, cdr3freqs, prodfreqs, reassignfreqs)
