from __future__ import print_function
import sys
import sqlite3
import random
from CreateAndImportClonesSqlite import *

if len(sys.argv) < 2:
    sys.exit("Usage: fix-multiple-V-assignments.py all_info.csv file(s)")

# Input
datafiles = sys.argv[1:]

def fixMultipleAssignedVGenes (datafile):
    # count entries
    nr_of_corrected_accs = 0
    nr_of_unique_accs = 0

    con = sqlite3.connect(":memory:")
    # con = sqlite3.connect("fix.db")
    cur = con.cursor()

    table = "all_info"
    basename = datafile.split("/")[-1]
    outfile = basename + "-clones-v-correction.txt"
    create_and_import(con, cur, table, datafile)

    # Get total number unique accession codes
    query = "select count(distinct acc) from all_info where V_sub!='None' AND J_sub!='None'"
    result = cur.execute(query)
    for row in result:
        nr_of_unique_accs = int(row[0])

    # Get number unique accession codes where V gene has to be re-assigned
    query = "select count(distinct acc) from all_info where cast(nr_v_subs as integer)>1"
    result = cur.execute(query)
    for row in result:
        nr_of_corrected_accs = int(row[0])

    # Stop if there are no accession codes with multiple assigned V genes
    if nr_of_corrected_accs == 0:
        return(nr_of_corrected_accs, nr_of_unique_accs)

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
        # print(query)
        result = cur.execute(query)
        cdr3peps = list()
        V_subs = list()
        freqs = list()
        for row in result:
            row = [str(c) for c in row]  # convert to string
            cdr3peps.append(row[0])
            V_subs.append(row[1])
            freqs.append(int(row[2]))

        # If there are no entries were unique V genes were assigned retrieve the info from accessions with multiple assigned V genes
        if len(freqs) == 0:
            query = "select cdr3pep,V_sub,count(distinct acc) from all_info where cdr3pep='" + cdr3pep + "' group by cdr3pep,V_sub order by count(distinct acc) desc"
            # print(query)
            result = cur.execute(query)
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

            # print(cdr3pep, acc, cdr3peptides[cdr3pep]['V_sub'])
            # print(V_subs)
            # print(fractions_cumulatief)
            # print("Assigned:", assign, assign_v)

            # Update the all_info table, set the newly assigned V gene for the accession code
            query = "UPDATE all_info SET V_sub='" + assign_v + "' WHERE acc='" + acc + "'"
            print(query)
            cur.execute(query)
    con.commit()

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
    fhClonesSubs = open(outfile, "w")
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

    # Close database connection
    con.close()

    return(nr_of_corrected_accs, nr_of_unique_accs)


########## MAIN ###############

datafiles = ["/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/PS043_S135_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv","/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/S074-129_S129_L001.assembled-CGATCGAT-IGH_HUMAN-all_info.csv"]

for datafile in datafiles:
    (nr_of_corrected_accs, nr_of_unique_accs) = fixMultipleAssignedVGenes(datafile)
    print("Corrected and total reads:", nr_of_corrected_accs, nr_of_unique_accs)
