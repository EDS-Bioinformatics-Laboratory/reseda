from __future__ import print_function
import sys
import sqlite3

# import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt

# Input
if len(sys.argv)<2:
    sys.exit("Usage: topX-reads-umis.py *-clones-subs.csv")
clonesFiles = sys.argv[1:]

########## Functions ###########

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
        cur.execute(query, row)
    con.commit()
    fh.close()

def compareTop (n):
    '''
    Description: get top N clones based on reads and on umi's. Check how much overlap there is between the lists
    In: n (get top N clones)
    Out: percentage overlap
    '''

    # Get top X based on read count
    top_based_on_reads = list()
    query = "select * from clones order by cast(freq as integer) desc limit " + str(n)       # top X reads
    result = cur.execute(query)
    for row in result:
        (v,j,cdr3) = (str(row[0]), str(row[1]), str(row[2]))
        top_based_on_reads.append("-".join([v,j,cdr3]))

    # Get top X based on UMI count
    top_based_on_umis = list()
    query = "select * from clones order by cast(uniq_umis as integer) desc limit " + str(n)  # top X umis
    result = cur.execute(query)
    for row in result:
        (v,j,cdr3) = (str(row[0]), str(row[1]), str(row[2]))
        top_based_on_umis.append("-".join([v,j,cdr3]))

    # Get percentage overlapping clones
    in_both = [val for val in top_based_on_umis if val in top_based_on_reads]
    perc_in_common = 100.0 * len(in_both) / n

    return(perc_in_common)

def makeBarChart (x,y):
    x_pos = np.arange(len(x))
    
    fig, ax = plt.subplots()
    ax.bar(x_pos, y, align='center', alpha=0.5)
    plt.xticks(x_pos, x)
    ax.set_xlabel('Top X clones')
    ax.set_ylabel('Common clones in topX UMIs versus reads (%)')
    ax.set_title(sampleName)
    
    try:
        fig.savefig(plotfile)
        print("Wrote", plotfile, "to disk")
    except:
        sys.exit("cannot write plotfile to disk")

########### Main ###########

for clonesFile in clonesFiles:
    # clonesFile = "/mnt/immunogenomics/RUNS/run04-20151116-miseq/test-reads-vs-umis/SP-CB16_S49_L001.assembled-ATGCATGC-IGH_HUMAN-clones-subs.csv"
    sampleName = clonesFile.split("/")[-1]
    plotfile = sampleName + "-topX.png"

    con = sqlite3.connect(":memory:")
    #con = sqlite3.connect("test.db")
    cur = con.cursor()

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

    # Determine top X clones based on read-frequency and on umi-count
    x = [1,5,10,25,50,100]
    y = list()
    for n in x:
        y.append(compareTop(n))

    # Make barchart
    makeBarChart(x,y)

    con.close()
