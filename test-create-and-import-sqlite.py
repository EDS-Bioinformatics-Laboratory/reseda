from __future__ import print_function
import sqlite3
from CreateAndImportClonesSqlite import *

def importData (datafile):
    # load all_info into database

    table = "aa_reads"
    create_and_import(con, cur, table, datafile)

#con = sqlite3.connect(":memory:")
con = sqlite3.connect("test.db")
cur = con.cursor()

datafile = "/home/narya/TMP/20110405_run81_paul_marieke/out_human_Bcells/aa/2016-04-15-AA.reads-81-BCRh-script20130529-.csv"
importData(datafile)

con.close()

print("Done")
