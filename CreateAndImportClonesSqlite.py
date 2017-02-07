from __future__ import print_function
# import sqlite3
import sys


def create_table(cur, name, colnames):
    '''
    Description: create a table in the database
    In: name (table name), colnames (list with column names)
    Out: query (string)
    '''

    query = "DROP TABLE IF EXISTS " + name + ";"
    print(query)
    cur.execute(query)

    query = "CREATE TABLE " + name + " ("
    query = query + ", ".join(colnames)
    query = query + ");"
    print(query)
    cur.execute(query)
    return(query)


def import_data(con, cur, datafile, delim, table, colnames):
    '''
    Description: Reads a csv file and inserts the content in a table
    In: datafile (path), delim (e.g. ',' or ' '), table (name of the table in database)
    Out: query (string)
    '''
    fh = open(datafile, 'r')
    fh.readline()  # skip first line
    for row in fh:
        row = row.rstrip()
        row = row.replace('"', '')
        row = row.split(delim)
        query = "INSERT INTO " + table + " ("
        query = query + ", ".join(colnames)
        query = query + ") VALUES ("
        query = query + ", ".join(len(colnames) * ["?"]) + ");"
        # print(query, row)
        cur.execute(query, row)
    con.commit()
    fh.close()


def create_and_import(con, cur, table, datafile):
    '''
    Description: create a new table and import the data
    In: str tablename, str filename
    '''
    delim = "\t"

    # Get header and close the file again
    try:
        fh = open(datafile, "rb")
    except:
        sys.exit("cannot open file" + datafile)
    header = fh.readline()
    header = header.rstrip()
    header = header.replace(":1", "2")
    header = header.replace(".", "_")
    header = header.replace('"', "")
    fh.close()

    # Create table and fill table with data
    colnames = header.split(delim)
    if colnames[0] == "":
        colnames[0] = "someId"
    if "group" in colnames:
        colnames[colnames.index("group")] = "grp"
    # colnames = ["blah"] + colnames   # only for AA_reads!!! Need to disable this line
    print(colnames)
    create_table(cur, table, colnames)
    import_data(con, cur, datafile, delim, table, colnames)

# ************ Main **************

# # Assume that script was imported when no arguments are given
# if len(sys.argv) == 1:
#     pass

# # Execute this when arguments are given to the script, else assume that the script is imported
# elif len(sys.argv) == 3:
#     #con = sqlite3.connect(":memory:")
#     con = sqlite3.connect("test.db")
#     cur = con.cursor()

#     #allInfoFile = "/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/PS043_S135_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv"
#     allInfoFile = sys.argv[1]
#     create_and_import(con, cur,"all_info", allInfoFile)

#     #clonesFile = "/mnt/immunogenomics/RUNS/run05-20151218-miseq/results-tbcell/final/correct-mid/PS043_S135_L001.assembled-ACTGACTG-TRB_HUMAN-clones-subs.csv"
#     clonesFile = sys.argv[2]
#     create_and_import(con, cur, "clones", clonesFile)

#     # Close the connection
#     con.close()

#     print("Finished")

#     exit()

# # Error message when there are arguments, but not the right amount
# else:
#     sys.exit("Usage: create_and_import.py path-to/all_info.csv path-to/clones-subs.csv")
