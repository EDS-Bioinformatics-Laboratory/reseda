from __future__ import print_function
import sqlite3
import random
from CreateAndImportClonesSqlite import *

#con = sqlite3.connect(":memory:")
con = sqlite3.connect("fix.db")
cur = con.cursor()

table = "all_info"
datafile = "/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/PS043_S135_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv"
#create_and_import(con, cur, table, datafile)

# Retrieve all cdr3peptides where multiple V's were assigned
query = "select cdr3pep,acc,V_sub from all_info where cast(nr_v_subs as integer)>1"
result = cur.execute(query)
cdr3peptides = dict()
for row in result:
    row = [str(c) for c in row]  # convert to string
    cdr3pep = row[0]
    acc = row[1]
    V_sub = row[2]

    cdr3peptides[cdr3pep] = cdr3peptides.get(cdr3pep, dict())  # make new sub dictionary if it doesn't exist
    cdr3peptides[cdr3pep]['acc'] = cdr3peptides[cdr3pep].get('acc', list()) + [acc]
    cdr3peptides[cdr3pep]['V_sub'] = cdr3peptides[cdr3pep].get('V_sub', list()) + [V_sub]


# Get clone frequency for every cdr3pep where the accessions got one unique V assigned
for cdr3pep in cdr3peptides:
    query = "select cdr3pep,V_sub,count(distinct acc) from all_info where cast(nr_v_subs as integer)=1 and cdr3pep='" + cdr3pep + "' group by cdr3pep,V_sub order by count(distinct acc) desc"
    result = cur.execute(query)
    cdr3peps = list()
    V_subs = list()
    freqs = list()
    for row in result:
        row = [str(c) for c in row]  # convert to string
        cdr3peps.append(row[0])
        V_subs.append(row[1])
        freqs.append(int(row[2]))

    # Calculate the fraction of clone frequency and the cumulative fraction
    total = float(sum(freqs))
    fractions = [fr/total for fr in freqs]
    fractions_cumulatief = [sum(fractions[:i+1]) for i in range(len(fractions))]
    if fractions_cumulatief[-1] != 1.0:
        fractions_cumulatief[-1] = 1.0

    for acc in cdr3peptides[cdr3pep]['acc']:
        # Generate random number between 0 and 1, assign a V gene based on fractions_cumulatief
        assign = random.random()
        for i in range(len(fractions_cumulatief)):
            if assign <= fractions_cumulatief[i]:
                assign_v = V_subs[i]
                break

        print(cdr3pep, acc, cdr3peptides[cdr3pep]['V_sub'])
        print(V_subs)
        print(fractions_cumulatief)
        print("Assigned:", assign, assign_v)

# Close the connection
con.close()
