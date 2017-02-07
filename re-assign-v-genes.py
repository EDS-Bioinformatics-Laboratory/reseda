from __future__ import print_function
import sys
import sqlite3
from CreateAndImportClonesSqlite import *
import matplotlib.pyplot as plt
import numpy as np


def importData(datafile, fhLog):
    # load all_info into database

    table = "all_info"
    # basename = datafile.split("/")[-1]
    create_and_import(con, cur, table, datafile)

    # Count entries with low quality CDR3
    result = cur.execute("SELECT COUNT(DISTINCT acc) FROM all_info WHERE CAST(cdr3_qual_min as int)<30")
    for row in result:
        low_qual_cdr3 = row[0]

    # Count entries with no V/J assigned
    result = cur.execute("SELECT COUNT(DISTINCT acc) FROM all_info WHERE V_gene='None' OR J_gene='None'")
    for row in result:
        no_vj = row[0]

    # Count entries with low quality CDR3 and no V/J assigned
    result = cur.execute("SELECT COUNT(DISTINCT acc) FROM all_info WHERE V_gene='None' OR J_gene='None' OR CAST(cdr3_qual_min as int)<30")
    for row in result:
        filtered_out = row[0]

    print(datafile, low_qual_cdr3, no_vj, filtered_out, file=fhLog)

    # # Remove all entries where the min_base_qual is below 30 or genes are not found
    # query = "DELETE FROM all_info WHERE V_gene='None' OR J_gene='None' OR CAST(cdr3_qual_min as int)<30"
    # print(query)
    # cur.execute(query)


def plotNrGenes(datafile, geneCol):
    basename = datafile.split("/")[-1]

    # query: get cdr3 and nr of different V or J genes
    query = "select cdr3pep,count(distinct " + geneCol + ") from all_info group by cdr3pep order by count(distinct " + geneCol + ") desc"
    result = cur.execute(query)
    data = list()
    for row in result:
        data.append(int(row[1]))

    if len(data) == 0:
        return

    # Start figure
    fig = plt.figure()
    fig.add_subplot(1, 1, 1)

    # Make histogram
    hist, bins = np.histogram(data, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.xticks(np.arange(min(data) - 2, max(data) + 1, 1))
    # ax.set_yscale('log')
    plt.bar(center, hist, align='center', width=width)
    plt.savefig(basename + ".rr.hist.svg")


def getCDR3WithMultipleGenes(geneCol):
    # query: get cdr3 and nr of different V or J genes with a count > 1
    query = "select cdr3pep,count(distinct " + geneCol + "),count(distinct acc) from all_info group by cdr3pep having count(distinct " + geneCol + ")>1 order by count(distinct " + geneCol + ") desc"
    result = cur.execute(query)
    cdr3s = list()
    for row in result:
        cdr3s.append(str(row[0]))
    return(cdr3s)


def getFreqGenes(cdr3, geneCol):

    # get the V or J genes
    query = "select " + geneCol + ",count(distinct acc) from all_info where cdr3pep='" + cdr3 + "' group by cdr3pep," + geneCol + " order by count(distinct acc) desc"
    result = cur.execute(query)
    gene = list()
    freq = list()
    for row in result:
        gene.append(str(row[0]))
        freq.append(int(row[1]))

    return(gene, freq)


def fraction(freqs):
    total = float(sum(freqs))
    frac = [f / total for f in freqs]
    frac_cum = [sum(frac[:i + 1]) for i in range(len(frac))]
    return(frac, frac_cum)


def getOffset(cutOff, frac_cum):
    return([e > cutOff for e in frac_cum].index(True))


# datafiles = ["/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/PS043_S135_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/S074-129_S129_L001.assembled-CGATCGAT-IGH_HUMAN-all_info.csv"]
# datafiles += ["/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB16_S49_L001.assembled-ATGCATGC-IGH_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB19_S52_L001.assembled-CGATCGAT-IGH_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB21_S54_L001.assembled-ACGTACGT-IGH_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/UNIFI68-9N-MID1-Exo_S162_L001.assembled-ACGTACGT-TRB_HUMAN-all_info.csv"]

# datafiles = ["/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/PS043_S135_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv"]
# datafiles = ["/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB16_S49_L001.assembled-ATGCATGC-IGH_HUMAN-all_info.csv"]
if __name__ == '__main__':

    # Input
    if len(sys.argv) < 2:
        sys.exit("Usage: re-assign-v-genes.py all_info.csv file(s)")
    datafiles = sys.argv[1:]

    cutOff = 0.7

    geneCol = "V_sub"

    fhOut = open("log-fix-multiple-" + geneCol + "-assignments.txt", "w")
    print("datafile corrected_accessions total_accessions", file=fhOut)
    for datafile in datafiles:
        outfile = datafile.split("/")[-1] + ".rr.csv"
        allfile = datafile.split("/")[-1] + ".rr.all_info.csv"
        clonefile = datafile.split("/")[-1] + ".rr.clones_subs.csv"

        try:
            fhOut = open(outfile, "w")
            print("cdr3\tnew\tgenes\tfreqs\tfrac", file=fhOut)
            fhAllInfo = open(allfile, "w")
            fhClonesSubs = open(clonefile, "w")
            fhLog = open(datafile.split("/")[-1] + ".quality-filter.log", "w")
            print("datafile low-cdr3 no-VJ filtered-out", file=fhLog)
        except:
            sys.exit("cannot write to disk")

        con = sqlite3.connect(":memory:")
        # con = sqlite3.connect("rr.db")
        cur = con.cursor()

        importData(datafile, fhLog)
        # plotNrGenes (datafile, geneCol)
        cdr3s = getCDR3WithMultipleGenes(geneCol)

        for cdr3 in cdr3s:
            # get V's, nr of accessions with that V (frequency)
            (genes, freqs) = getFreqGenes(cdr3, geneCol)
            # calculate fraction and cumulative fraction of frequency
            (frac, frac_cum) = fraction(freqs)
            # determine which V's to combine (if the first V occurs more than 70% take that one, if first two V's occur more than 70% combine these, etc)
            i = getOffset(cutOff, frac_cum)
            new = genes[0:i + 1]
            new.sort()
            new = "+".join(new)
            genes = ",".join(genes)
            freqs = ",".join([str(c) for c in freqs])
            frac = ",".join([str(c) for c in frac])
            print("\t".join([cdr3, new, genes, freqs, frac]), file=fhOut)
            # update all_info table
            query = "update all_info set " + geneCol + "='" + new + "' where cdr3pep='" + cdr3 + "'"
            print(query)
            cur.execute(query)

        # Write all_info to a file
        fhAllInfo = open(allfile, "w")
        result = cur.execute('SELECT * FROM all_info')
        print("\t".join([description[0] for description in result.description]), file=fhAllInfo)
        for row in result:
            str_row = [str(c) for c in row]
            print("\t".join(str_row), file=fhAllInfo)
        fhAllInfo.close()

        # Create a clone report based on V_sub, J_sub and CDR3peptide
        query = "DROP TABLE IF EXISTS clones_subs"
        print(query)
        cur.execute(query)
        query = "CREATE TABLE clones_subs AS SELECT V_sub, J_sub, cdr3pep, count(DISTINCT acc) AS freq, count(DISTINCT beforeMID) AS uniq_umis FROM all_info WHERE V_sub!='None' AND J_sub!='None' GROUP BY V_sub, J_sub, cdr3pep"
        print(query)
        cur.execute(query)

        result = cur.execute("SELECT SUM(freq) AS total_reads, SUM(uniq_umis) AS total_umis FROM clones_subs")
        for row in result:
            try:
                total_reads = int(row[0])
            except:
                total_reads = 0
            try:
                total_umis = int(row[1])
            except:
                total_umis = 0

        # Write clones to a file
        result = cur.execute('SELECT * FROM clones_subs ORDER BY freq DESC')
        print("\t".join([description[0] for description in result.description]) + "\tread_perc\tumi_perc", file=fhClonesSubs)
        for row in result:
            str_row = list()
            for i in range(len(row)):
                str_row.append(str(row[i]))
            # Calculate percentage
            read_perc = 100 * float(str_row[3]) / float(total_reads)  # column 3 is the frequency
            str_row.append(str(read_perc))
            umi_perc = 100 * float(str_row[4]) / float(total_umis)  # column 3 is the frequency
            str_row.append(str(umi_perc))
            print("\t".join(str_row), file=fhClonesSubs)

        con.close()
        fhOut.close()
        fhAllInfo.close()
        fhClonesSubs.close()
        fhLog.close()
