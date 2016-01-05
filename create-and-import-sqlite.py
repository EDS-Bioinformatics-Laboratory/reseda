from __future__ import print_function
import sqlite3
import sys


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
    header = fh.readline()  # skip first line
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
    fh.close()


########### Main ###########

#con = sqlite3.connect(":memory:")
con = sqlite3.connect("test.db")
cur = con.cursor()

allInfoFile = "/mnt/immunogenomics/RUNS/run04-20151116-miseq/test-reads-vs-umis/SP-CB16_S49_L001.assembled-ATGCATGC-IGH_HUMAN-all_info.csv"

# Get header and close the file again
try:
    fh = open(allInfoFile, "rb")
except:
    sys.exit("cannot open file")
header = fh.readline()
header = header.replace(":1", "2")
fh.close()

# Create table and fill table with data
colnames = header.split()
create_table("all_info",colnames)
import_data(allInfoFile, "\t", "all_info", colnames)

clonesFile = "/mnt/immunogenomics/RUNS/run04-20151116-miseq/test-reads-vs-umis/SP-CB16_S49_L001.assembled-ATGCATGC-IGH_HUMAN-clones-subs.csv"

# Get header and close the file again
try:
    fh = open(clonesFile, "rb")
except:
    sys.exit("cannot open file")
header = fh.readline()
fh.close()

# Create table and fill table with data
colnames = header.split()
create_table("clones",colnames)
import_data(clonesFile, "\t", "clones", colnames)

# Close the connection
con.close()

print("Finished")

exit()
