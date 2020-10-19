
# coding: utf-8

# In[1]:


import pandas as pd
import os


# In[2]:


# Input: CDR3 files of all Roche runs (except run214 and run258)
mydir = "/mnt/immunogenomics/RUNS/runXXX-roche-LN2/CDR3_20201015/"
myfiles = [x for x in os.listdir(mydir) if x.endswith("-IGH_HUMAN-CDR3.csv")]
myfiles[:10]


# In[3]:


# Input: list of CDR3's that I received from Aram
df_cdr3 = pd.read_csv("2020-08-21.presyn-top100-clones-contamination-table-on-pt.csv", sep=",")
df_cdr3.columns = ['unnamed', 'CDR3_x', 'Samples']
df_cdr3.head()


# In[4]:


# Output
# Files will be written in the current directory with this name: outfile = myfile.replace(".csv", "-matches.csv")


# In[5]:


f = lambda x: x.replace(" x", "")
df_cdr3['CDR3'] = [x for x in map(f, df_cdr3['CDR3_x'])]
df_cdr3.head()


# ## Find the CDR3s in all Roche runs

# In[6]:


for myfile in myfiles:
    try:
        df = pd.read_csv(mydir + myfile, sep="\t", header=None)
        df.columns = ['acc', 'frame', 'CDR3pep', 'CDR3nuc', 'CDR3pep-short', 'CDR3nuc-short', 'qc1', 'qc2', 'qc3', 'basecalls', 'x', 'y', 'z']
    except:
        print("ERROR: couldn't read", mydir + myfile)
        continue
    df = pd.merge(df, df_cdr3, how="inner", left_on="CDR3pep", right_on="CDR3")
    outfile = myfile.replace(".csv", "-matches.csv")
    df.to_csv(outfile)
    print("Wrote", outfile, "to disk")

