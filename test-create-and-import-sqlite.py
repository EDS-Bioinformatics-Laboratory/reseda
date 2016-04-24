from __future__ import print_function
import sqlite3
from CreateAndImportClonesSqlite import *

def importData (tablename, datafile):
    # load all_info into database

    create_and_import(con, cur, tablename, datafile)

#con = sqlite3.connect(":memory:")
con = sqlite3.connect("run06.db")
cur = con.cursor()

# importData("aa_reads", "/home/narya/TMP/20110405_run81_paul_marieke/out_human_Bcells/aa/2016-04-15-AA.reads-81-BCRh-script20130529-.csv")
importData("clones", "/home/barbera/TMP/run06-clones_subs.csv")
con.close()

print("Done")
