
# coding: utf-8

# In[1]:


import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


# In[2]:


option_comrev = 1  # 0: do not make complent reverse, 1: complement reverse the experiment sequences
option_project_name = "DO-051V1"
option_outdir = "./"
# option_vgene = "top"
# option_jgene = "top"
option_vgene = "IGHV3-23*01"
option_jgene = "IGHJ4*02"


# In[3]:


mydir = "20210909-run42-DOMINO-IGH/"
myfiles = sorted([x for x in os.listdir(mydir) if x.startswith("DO-051") and x.endswith("all_info.csv")])
myfiles


# In[4]:


def readAllInfo(mydir, f):
    samplename, rest = f.split("_L001")
    df = pd.read_csv(mydir + f, sep="\t")
    df["Sample"] = samplename
    return(df)
df = readAllInfo(mydir, myfiles[0])
for f in myfiles[1:]:
    df = pd.concat([df, readAllInfo(mydir, f)])
df.head()


# In[5]:


df.columns


# In[6]:


# Which VJ combination occurs most frequent?
df_vj_count = df.groupby(['V_gene', 'J_gene']).agg({'Sample': 'nunique', 'acc': 'nunique'})
df_vj_count = df_vj_count.sort_values(by=['Sample','acc'], ascending=False)
df_vj_count = df_vj_count.reset_index()
df_vj_count.head()


# In[7]:


if option_vgene != "top" and option_jgene: # Were both the V and J genes given as an option? Select those
    v_top = option_vgene
    j_top = option_jgene
else:                                      # Else, select the VJ combination that occurs often
    v_top = df_vj_count['V_gene'][0]
    j_top = df_vj_count['J_gene'][0]


# In[8]:


# Retrieve sequences for VJ combination and make these sequences unique for further analysis
concatenate = lambda x: "|".join(list(set(x)))
df_selection = df[(df["V_gene"] == v_top) & (df["J_gene"] == j_top)]
df_selection = df_selection.groupby("seq").agg({'acc': concatenate, 'Sample': concatenate}).reset_index()
print("entries:", len(df_selection))
df_selection.head()


# In[9]:


# Open files for writing the fasta sequences
out_fastaV = option_outdir + option_project_name + ".V.fasta"
out_fastaJ = option_outdir + option_project_name + ".J.fasta"

fhOutV = open(out_fastaV, "w")
fhOutJ = open(out_fastaJ, "w")


# In[10]:


# Retrieve the reference V and J sequences
v_seq, j_seq = "", ""
for record in SeqIO.parse(open("../reference/IGHV_human.fasta"), "fasta"):
    if v_top in record.id:
        v_seq = str(record.seq).upper()
for record in SeqIO.parse(open("../reference/IGHJ_human.fasta"), "fasta"):
    if j_top in record.id:
        j_seq = str(record.seq).upper()

print(">" + v_top)
print(v_seq)
print(">" + j_top)
print(j_seq)
print(">" + v_top, file=fhOutV)
print(v_seq, file=fhOutV)
print(">" + j_top, file=fhOutJ)
print(j_seq, file=fhOutJ)


# In[11]:


comrev = lambda s: str(Seq(s).reverse_complement()).upper()
comrev("AAAACGATCGATCGTATCGCCTCCCTCGCGCCATGTGTTTCCC")


# In[12]:


for i in df_selection.index:
    print(">" + df_selection.iloc[i]['acc'] + "|" + df_selection.iloc[i]['Sample'], file=fhOutV)
    if option_comrev == 1:
        print(comrev(df_selection.iloc[i]['seq']), file=fhOutV)
    else:
        print(df_selection.iloc[i]['seq'], file=fhOutV)
        
    print(">" + df_selection.iloc[i]['acc'], file=fhOutJ)
    if option_comrev == 1:
        print(comrev(df_selection.iloc[i]['seq']), file=fhOutJ)
    else:
        print(df_selection.iloc[i]['seq'], file=fhOutJ)


# In[13]:


fhOutV.close()
fhOutJ.close()
print("Wrote", out_fastaV, "to disk")
print("Wrote", out_fastaJ, "to disk")

