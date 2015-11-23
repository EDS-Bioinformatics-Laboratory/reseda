from __future__ import print_function
import sys

infile = "/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/S074-135_S12_L001.assembled-ACTGACTG-IGH_HUMAN-all_info.csv"

def analyzeCdr3Qual (infile):
    filepath = infile.split("/")
    filename = filepath[-1]
    outMin = filename + ".min_qual.csv"
    outAvg = filename + ".avg_qual.csv"

    try:
        fh = open(infile, "rb")
    except:
        sys.exit("cannot open file")

    firstLine = fh.readline()
    colnames = firstLine.split()
    i_cdr3_qual_min = 7  # default col name
    i_cdr3_qual_avg = 9
    for i in range(len(colnames)):
        if colnames[i] == "cdr3_qual_min":
            i_cdr3_qual_min = i
        elif colnames[i] == "cdr3_qual_avg":
            i_cdr3_qual_avg = i

    min_qual = dict()
    avg_qual = dict()
    total_reads = 0
    discarded_reads = 0
    for line in fh:
        line = line.rstrip()
        col = line.split()

        cdr3_qual_min = int(float(col[i_cdr3_qual_min]))
        cdr3_qual_avg = int(float(col[i_cdr3_qual_avg]))

        min_qual[cdr3_qual_min] = min_qual.get(cdr3_qual_min, 0) + 1
        avg_qual[cdr3_qual_avg] = avg_qual.get(cdr3_qual_avg, 0) + 1

        total_reads = total_reads + 1
        if cdr3_qual_min <= 15:
            discarded_reads = discarded_reads + 1

    fh.close()

    try:
        fhMin = open(outMin, "w")
        fhAvg = open(outAvg, "w")
    except:
        sys.exit("cannot write files")

    for key in sorted(min_qual, reverse=True):
        print(key, min_qual[key], file=fhMin)

    for key in sorted(avg_qual, reverse=True):
        print(key, avg_qual[key], file=fhAvg)

    perc_discarded = 100.0 * discarded_reads / total_reads
    print("CDR3 with min_qual below/equal to 15:", discarded_reads, "/", total_reads, "=", perc_discarded, "%")

    fhMin.close()
    fhAvg.close()

############ Main #########

analyzeCdr3Qual(infile)
