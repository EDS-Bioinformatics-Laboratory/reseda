from __future__ import print_function
import sys
import sqlite3

# import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.pyplot as plt

########## Functions ###########

def getSamplesOverlappingHECs (threshold):
    '''
    Description: get HECs based on read_perc threshold, count how many samples show overlap
    In: threshold (float or int)
    Out: number of samples with overlapping clones
    '''

    query = "drop table if exists tmp"
    result = cur.execute(query)

    query = "create table tmp as select V_sub, J_sub, cdr3pep, count(distinct Sample) as nr_samples from clones where cast(read_perc as numeric)>"+ str(threshold) +" group by V_sub, J_sub, cdr3pep having nr_samples>1 order by nr_samples desc"
    result = cur.execute(query)

    query = "drop table if exists tmp_clones"
    result = cur.execute(query)

    query = "create table tmp_clones as select b.* from tmp as a,clones as b where a.V_sub=b.V_sub and a.J_sub=b.J_sub and a.cdr3pep=b.cdr3pep and cast(b.read_perc as numeric)>" + str(threshold)
    result = cur.execute(query)

    return(1)

    # # Get top X based on UMI count
    # top_based_on_umis = list()
    # query = "select * from clones order by cast(uniq_umis as integer) desc limit " + str(n)  # top X umis
    # result = cur.execute(query)
    # for row in result:
    #     (v,j,cdr3) = (str(row[0]), str(row[1]), str(row[2]))
    #     top_based_on_umis.append("-".join([v,j,cdr3]))

    # return(perc_in_common)

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

# con = sqlite3.connect(":memory:")
con = sqlite3.connect("run06.db")
cur = con.cursor()

# Determine top X clones based on read-frequency and on umi-count
# x = [0.1,0.25,0.5,0.75,1,2,5]
x = [0.5]
y = list()
for t in x:
    y.append(getSamplesOverlappingHECs(t))

# Make barchart
# makeBarChart(x,y)

con.close()
