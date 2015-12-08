from __future__ import print_function
import sqlite3
import sys

# if len(sys.argv) < 8:
#     sys.exit("Usage: combine-immuno-data.py midFile cdr3File vFile jFile seqFile outFile clonesFile clonesGroupedFile totalFile")

# # Input files
# [midFile,cdr3File,vFile,jFile,seqFile,outFile,cloneFile,cloneGroupedFile,totalFile] = sys.argv[1:10]

# Input files - TEST
midFile = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-tbcell/reports/BCRh_S40_L001.assembled-report.txt"
cdr3File = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/paul/BCRh_S40_L001.assembled.fastq.gz-IGH_HUMAN-CDR3.csv"
vFile = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/BCRh/BCRh_S40_L001.assembled-IGHV_human-easy-import.txt"
jFile = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/BCRh/BCRh_S40_L001.assembled-IGHJ_human-easy-import.txt"
seqFile = "/mnt/immunogenomics/RUNS/run03-20150814-miseq/results-pear/paul/BCRh_S40_L001.assembled.fastq.gz-IGH_HUMAN.csv"
outFile = "all_info.txt"
cloneFile = "clones.txt"
cloneGroupedFile = "clones-grouped.txt"
totalFile = "total.txt"

# Output file
fhOut = open(outFile, 'w')
fhClones = open(cloneFile, 'w')
fhClonesGrouped = open(cloneGroupedFile, 'w')
fhTotal = open(totalFile, 'w')

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

def clean_name (gene):
    '''
    Description: cleans up unnessary info in the gene name that was returned by the aligner, also removes allele info
    In: V or J name (e.g. 'M99672|IGHV3-43*01|Homo')
    Out: cleaned up gene name (e.g. 'IGHV3-43')
    '''
    tmp = gene.split("|")

    if len(tmp) == 3:
        gene = tmp[1]   # remove accession at start and 'Homo' or 'Mouse' at the end
        gene = gene.split("*")[0]   # allele info
    else:
        return(gene)


    return(gene)

########### Main ###########

# con = sqlite3.connect(":memory:")
con = sqlite3.connect("test.db")
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
colnames = ["acc","readingframe_seq","seq","pep","qual"]
create_table("seq",colnames)
import_data(seqFile, "\t", "seq", colnames)

# Combine all tables into one big table: all_info
query = "CREATE TABLE all_info AS SELECT * FROM mid JOIN cdr3 USING (acc) LEFT OUTER JOIN v USING (acc) LEFT OUTER JOIN j USING (acc) LEFT OUTER JOIN seq USING (acc);"
print(query)
cur.execute(query)

# Dump info from this table to a file
result = cur.execute('SELECT * FROM all_info')
header = "\t".join([description[0] for description in result.description])
#print(header + "\tV\tJ", file=fhOut)
colnames = header.split()
i_v = 0
i_j = 0
# Check which column contains the V_gene
for i in range(len(colnames)):
    if colnames[i] == "V_gene":
        i_v = i
    elif colnames[i] == "J_gene":
        i_j = i
# Print table to file
for row in result:
    str_row = list()
    for i in range(len(row)):
        str_row.append(str(row[i]))

    # Make changes to the V_gene name
    V = clean_name(str_row[i_v])
    V_main = V.split("-")[0]
    J = clean_name(str_row[i_j])

    # Print the entry
    print("\t".join(str_row) + "\t" + V + "\t" + J + "\t" + V_main, file=fhOut)
fhOut.close()

# Import table all_info again and then create clone list
colnames.append("V")
colnames.append("J")
colnames.append("V_main")
query = "DELETE FROM all_info"
print(query)
cur.execute(query)
query = "ALTER TABLE all_info ADD COLUMN V"
print(query)
cur.execute(query)
query = "ALTER TABLE all_info ADD COLUMN J"
print(query)
cur.execute(query)
query = "ALTER TABLE all_info ADD COLUMN V_main"
print(query)
cur.execute(query)
import_data(outFile, "\t", "all_info", colnames)

# Create a clone report based on V, J and CDR3peptide
query = "CREATE TABLE clones AS SELECT V, J, cdr3pep, count(DISTINCT acc) AS freq FROM all_info WHERE V!='None' AND J!='None' GROUP BY V, J, cdr3pep"
print(query)
cur.execute(query)

# Write clones to a file
result = cur.execute('SELECT * FROM clones ORDER BY freq DESC')
print("\t".join([description[0] for description in result.description]), file=fhClones)
for row in result:
    str_row = list()
    for i in range(len(row)):
        str_row.append(str(row[i]))
    print("\t".join(str_row), file=fhClones)
fhClones.close()

# Create a clone report based on V_main, J and CDR3peptide
query = "CREATE TABLE clones_grouped AS SELECT V_main, J, cdr3pep, count(DISTINCT acc) AS freq FROM all_info WHERE V_main!='None' AND J!='None' GROUP BY V_main, J, cdr3pep"
print(query)
cur.execute(query)

# Write clones to a file
result = cur.execute('SELECT * FROM clones_grouped ORDER BY freq DESC')
print("\t".join([description[0] for description in result.description]), file=fhClonesGrouped)
for row in result:
    str_row = list()
    for i in range(len(row)):
        str_row.append(str(row[i]))
    print("\t".join(str_row), file=fhClonesGrouped)
fhClonesGrouped.close()

### Totals ###

result = cur.execute('select count(*) from all_info')
for row in result:
    print("Total rows in all_info:", row[0], file=fhTotal)

result = cur.execute('select count(distinct acc) from all_info')
for row in result:
    print("Unique reads in all_info:", row[0], file=fhTotal)

result = cur.execute("select count(*) from all_info where V_gene!='None' and J_gene!='None'")
for row in result:
    print("Total rows with V and J in all_info:", row[0], file=fhTotal)

result = cur.execute("select count(distinct acc) from all_info where V_gene!='None' and J_gene!='None'")
for row in result:
    print("Unique reads with V and J in all_info:", row[0], file=fhTotal)

result = cur.execute('SELECT SUM(freq) FROM clones')
for row in result:
    print("Total reads in clones table:", row[0], file=fhTotal)

fhTotal.close()

# Close the connection
con.close()

print("Finished")

exit()
