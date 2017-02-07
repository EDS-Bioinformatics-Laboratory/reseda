from __future__ import print_function
import sqlite3
from CreateAndImportClonesSqlite import *


def importData(tablename, datafile):
    # load all_info into database

    create_and_import(con, cur, tablename, datafile)


con = sqlite3.connect(":memory:")
# con = sqlite3.connect("run03.db")
cur = con.cursor()

# importData("aa_reads", "/home/narya/TMP/20110405_run81_paul_marieke/out_human_Bcells/aa/2016-04-15-AA.reads-81-BCRh-script20130529-.csv")
# importData("clones", "/home/barbera/TMP/run06-clones_subs.csv")
importData("all_info", "/home/barbera/TMP/HLA-BCRh_S36_L001.assembled-AGCTAGCT-IGH_HUMAN-all_info.csv.rr.all_info.csv.IGHV3-7-IGHJ6.csv")

# Create a clone report based on V_sub, J_sub and CDR3peptide
fhClonesSubs = open("clones.csv", "w")

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
fhClonesSubs.close()

con.close()
print("Done")
