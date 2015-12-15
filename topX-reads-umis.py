from __future__ import print_function
import sqlite3
import sys


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

########### Main ###########

#con = sqlite3.connect(":memory:")
con = sqlite3.connect("test.db")
cur = con.cursor()

for n in [1, 5, 10, 25, 50, 100]:
    print(n, compareTop(n))
