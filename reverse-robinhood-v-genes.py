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


datafiles = ["/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/PS043_S135_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/S074-129_S129_L001.assembled-CGATCGAT-IGH_HUMAN-all_info.csv"]
datafiles += ["/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB16_S49_L001.assembled-ATGCATGC-IGH_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB19_S52_L001.assembled-CGATCGAT-IGH_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB21_S54_L001.assembled-ACGTACGT-IGH_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/UNIFI68-9N-MID1-Exo_S162_L001.assembled-ACGTACGT-TRB_HUMAN-all_info.csv"]

# datafiles = ["/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/PS043_S135_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv"]
# datafiles = ["/mnt/immunogenomics/RUNS/run04-20151116-miseq/results-tbcell/final/correct-mid/SP-CB16_S49_L001.assembled-ATGCATGC-IGH_HUMAN-all_info.csv"]

fhOut = open("log-fix-multiple-V-assignments.txt", "w")
print("datafile corrected_accessions total_accessions", file=fhOut)
for datafile in datafiles:
    outfile = datafile.split("/")[-1] + ".rr.csv"
    allfile = datafile.split("/")[-1] + ".rr.all_info.csv"
    clonefile = datafile.split("/")[-1] + ".rr.clones_subs.csv"
    try:
        fhOut = open(outfile, "w")
        print("cdr3\tnew_v\tv_genes\tfreqs\tfrac", file=fhOut)
        fhAllInfo = open(allfile, "w")
        fhClonesSubs = open(clonefile, "w")
    except:
        sys.exit("cannot write to disk")

    con = sqlite3.connect(":memory:")
    # con = sqlite3.connect("rr.db")
    cur = con.cursor()

    importData(datafile)
    cdr3s = getCdr3WithMultipleVgenes()

    for cdr3 in cdr3s:
        # get V's, nr of accessions with that V (frequency)
        (v_genes, freqs) = getFreqVgenes(cdr3)
        # calculate fraction and cumulative fraction of frequency
        (frac, frac_cum) = fraction(freqs)
        # determine which V's to combine (if the first V occurs more than 70% take that one, if first two V's occur more than 70% combine these, etc)
        i = getOffset(0.7, frac_cum)
        new_v = v_genes[0:i+1]
        new_v.sort()
        new_v = "+".join(new_v)
        v_genes = ",".join(v_genes)
        freqs = ",".join([str(c) for c in freqs])
        frac = ",".join([str(c) for c in frac])
        print("\t".join([cdr3, new_v, v_genes, freqs, frac]), file=fhOut)
        # update all_info table
        query = "update all_info set V_sub='" + new_v + "' where cdr3pep='" + cdr3 + "'"
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
        total_reads = int(row[0])
        total_umis = int(row[1])

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
