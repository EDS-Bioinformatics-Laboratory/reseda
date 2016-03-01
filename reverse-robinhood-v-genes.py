from __future__ import print_function
import sys
import sqlite3
from CreateAndImportClonesSqlite import *

# # Input
# if len(sys.argv) < 2:
#     sys.exit("Usage: fix-multiple-V-assignments.py all_info.csv file(s)")
# datafiles = sys.argv[1:]

def importData (datafile):
    # load all_info into database

    table = "all_info"
    basename = datafile.split("/")[-1]
    create_and_import(con, cur, table, datafile)

def getCdr3WithMultipleVgenes ():
    # query: get cdr3 and nr of different V genes
    query = "select cdr3pep,count(distinct V_sub),count(distinct acc) from all_info group by cdr3pep having count(distinct V_sub)>1 order by count(distinct V_sub) desc"
    result = cur.execute(query)
    cdr3s = list()
    for row in result:
        cdr3s.append(str(row[0]))
    return(cdr3s)

def getFreqVgenes (cdr3):

    # get the V genes
    query = "select V_sub,count(distinct acc) from all_info where cdr3pep='" + cdr3 + "' group by cdr3pep,V_sub order by count(distinct acc) desc"
    result = cur.execute(query)
    v_gene = list()
    freq = list()
    for row in result:
        v_gene.append(str(row[0]))
        freq.append(int(row[1]))

    return(v_gene,freq)

def fraction (freqs):
    total = float(sum(freqs))
    frac = [f/total for f in freqs]
    frac_cum = [sum(frac[:i+1]) for i in range(len(frac))]
    return(frac, frac_cum)

def getOffset (cutOff, frac_cum):
    return([e>cutOff for e in frac_cum].index(True))

####### MAIN ########


# datafiles = ["/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/PS043_S135_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/S074-129_S129_L001.assembled-CGATCGAT-IGH_HUMAN-all_info.csv"]
# datafiles += ["/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB16_S49_L001.assembled-ATGCATGC-IGH_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB19_S52_L001.assembled-CGATCGAT-IGH_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB21_S54_L001.assembled-ACGTACGT-IGH_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/UNIFI68-9N-MID1-Exo_S162_L001.assembled-ACGTACGT-TRB_HUMAN-all_info.csv"]
datafiles = ["/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/PS043_S135_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv"]

fhOut = open("log-fix-multiple-V-assignments.txt", "w")
print("datafile corrected_accessions total_accessions", file=fhOut)
for datafile in datafiles:
    # con = sqlite3.connect(":memory:")
    con = sqlite3.connect("rr.db")
    cur = con.cursor()

    # importData(datafile)
    cdr3s = getCdr3WithMultipleVgenes()

    for cdr3 in cdr3s:
        (v_genes, freqs) = getFreqVgenes(cdr3)
        (frac, frac_cum) = fraction(freqs)
        # print(cdr3, v_genes, freqs, frac, frac_cum)
        i = getOffset(0.7, frac_cum)
        print(cdr3, v_genes[0:i+1],frac[0:i+1])


    # determine which V gene occurs most
    # assign V gene with highest occurrence if above XX% (70?)
    #assign multiple V genes when it is difficult to assign (e.g. 50/50%)

    con.close()
