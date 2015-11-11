from __future__ import print_function
import sqlite3
import sys

if len(sys.argv) < 7:
    sys.exit("Usage: combine-immuno-data.py midFile cdr3File vFile jFile seqFile outFile")

# Input files
[midFile,cdr3File,vFile,jFile,seqFile,outFile] = sys.argv[1:7]

# Input files - TEST
# midFile = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/data/midsort/BCRh_S40_L001_R2_001-report.txt"
# cdr3File = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/paul/BCRh_S40_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv"
# vFile = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/BCRh/BCRh_S40_L001.assembled-IGHV_human-easy-import.txt"
# jFile = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/BCRh/BCRh_S40_L001.assembled-IGHJ_human-easy-import.txt"
# seqFile = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/paul/BCRh_S40_L001.assembled.fastq.gz-IGH_HUMAN.csv"
# outFile = "all_info.txt"

# Output file
fhOut = open(outFile, 'w')

######### Functions ########

def create_table (name, colnames):
    '''
    Description: create a table in the database
    In: name (table name), colnames (list with column names)
    Out: query (string)
    '''
    query = "CREATE TABLE " + name + " ("
    query = query + ", ".join(colnames)
    query = query + ");"
    print(query)
    cur.execute(query)
    return(query)

def import_data (datafile, delim, table, colnames):
    '''
    Description: Reads a csv file and inserts the content in a table
    In: datafile (path), delim (e.g. ',' or ' '), table (name of the table in database)
    Out: query (string)
    '''
    fh = open(datafile, 'r')
    for row in fh:
        row = row.rstrip()
        row = row.split(delim)
        query = "INSERT INTO " + table + " ("
        query = query + ", ".join(colnames)
        query = query + ") VALUES ("
        query = query + ", ".join(len(colnames) * ["?"]) + ");"
        # print(query, row)
        cur.execute(query, row)
    con.commit()


########### Main ###########

con = sqlite3.connect(":memory:")
#con = sqlite3.connect("test.db")
cur = con.cursor()

# MID
colnames = ["acc","beforeMID","MID","afterMID"]
create_table("mid",colnames)
import_data(midFile, " ", "mid", colnames)

# CDR3
colnames = ["acc","readingframe","cdr3pep","cdr3nuc", "cdr3_qual_min", "cdr3_qual_max", "cdr3_qual_avg","cdr3_qual"]
create_table("cdr3",colnames)
import_data(cdr3File, "\t", "cdr3", colnames)

# V genes
colnames = ["acc","V_flag","V_gene"]
create_table("v",colnames)
import_data(vFile, "\t", "v", colnames)

# J genes
colnames = ["acc","J_flag","J_gene"]
create_table("j",colnames)
import_data(jFile, "\t", "j", colnames)

# Sequences
colnames = ["acc","readingframe","seq","pep","qual"]
create_table("seq",colnames)
import_data(seqFile, "\t", "seq", colnames)

# Combine all tables into one big table: all_info
query = "CREATE TABLE all_info AS SELECT * FROM mid JOIN cdr3 USING (acc) LEFT OUTER JOIN v USING (acc) LEFT OUTER JOIN j USING (acc) LEFT OUTER JOIN seq USING (acc);"
print(query)
cur.execute(query)

# Dump info from this table to a file
result = cur.execute('SELECT * FROM all_info')
print("\t".join([description[0] for description in result.description]), file=fhOut)
for row in result:
    str_row = list()
    for i in range(len(row)):
        str_row.append(str(row[i]))
    print("\t".join(str_row), file=fhOut)

# Close the connection
con.close()

print("Finished")

exit()
