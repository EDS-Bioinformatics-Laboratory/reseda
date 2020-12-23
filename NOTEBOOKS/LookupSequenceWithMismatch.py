
# coding: utf-8

# In[18]:


import os
import pandas as pd
import regex


# In[2]:


# output files
outHitsWithMismatch = "basta-allinfo-mismatch.xlsx"


# ## Read auto-reactive sequences

# In[3]:


myfile = "autoantibody_sequences.xlsx"


# In[4]:


xl = pd.ExcelFile(myfile)
sheet_names = xl.sheet_names  # see all sheet names
print(sheet_names)


# In[5]:


def read_excel_sheet(myfile, sheetname):
    df = pd.read_excel(myfile, sheet_name=sheetname, header=None)
    df.columns = ['sequence']
    df['sheet'] = sheetname
    return(df)


# In[6]:


df_cdr3 = read_excel_sheet(myfile, sheet_names[0])
for sheetname in sheet_names[1:]:
    df_cdr3 = pd.concat([df_cdr3, read_excel_sheet(myfile, sheetname)])
df_cdr3 = df_cdr3.reset_index()
df_cdr3.tail()


# In[7]:


df_cdr3.head()


# In[8]:


# Remove the ">>some-name" rows
patternDel = ">>"
filter = df_cdr3['sequence'].str.contains(patternDel)
df_cdr3 = df_cdr3[~filter]
df_cdr3.head()


# In[9]:


def f1(x): # remove first two amino acids if it starts with YY
    if x.startswith("YY"):
        return(x[2:])
    else:
        return(x)

def f2(x): # remove everything after the "VT" if the pattern "VTVS" is in the sequence
    x = x.split("VTVS")
    if len(x) > 2:
        print("WARNING: multiple times VTVS in sequence", x)
    return(x[0] + "VT")
    
df_cdr3['sequence_without_YY'] = [x for x in map(f1, df_cdr3['sequence'])]
df_cdr3['cdr3pep'] = [x for x in map(f2, df_cdr3['sequence_without_YY'])]
df_cdr3.head()


# ## Read all info files and lookup sequences, allow for mismatches

# In[10]:


allinfo_files = [x for x in os.listdir(".") if x.endswith(".all_info.csv")]
allinfo_files[:10]


# In[16]:


cdr3_list = list(set(df_cdr3['sequence']))
print("number of sequences to lookup:", len(cdr3_list))


# In[84]:


motif = "(" + "|".join(cdr3_list) + "){e<=1}"
motif = motif.replace(" ", "") # remove whitespace
motif = motif.replace("||", "|")
print(motif)
p = regex.compile(motif, regex.BESTMATCH)


# In[96]:


def lookupSubSequence(pep):
    # motif and p are global variables
    m = p.search(str(pep))
    if m is not None:   # a match is found
        match = m.group()
        
        # Lookup the original sequence
        motif_reverse = "(" + match + "){e<=1}"
        p_reverse = regex.compile(motif_reverse, regex.BESTMATCH)
        m_reverse = p_reverse.search(str(motif))
        orig_pattern = m_reverse.group()
        return(orig_pattern + "," + match)
    else:
        return(None)


# In[97]:


def lookupSubSequencePerFile(df_allinfo):
    df_allinfo['hit'] = [x for x in map(lookupSubSequence, df_allinfo['pep'])]
    df_tmp = df_allinfo[df_allinfo['hit'].notna()]
    df_tmp["orig_seq"] = [x.split(",")[0] for x in df_tmp["hit"]]
    df_tmp["found_seq"] = [x.split(",")[1] for x in df_tmp["hit"]]
    return(df_tmp)


# In[98]:


df_allinfo = pd.read_csv(allinfo_files[0], sep="\t")
sample_name, rest = allinfo_files[0].split("_L001")
df_allinfo["Sample"] = sample_name
df_lookup_long = lookupSubSequencePerFile(df_allinfo)
print(sample_name, len(df_lookup_long))


# In[ ]:


for allinfo_file in allinfo_files[1:]:
    df_allinfo = pd.read_csv(allinfo_file, sep="\t")
    sample_name, rest = allinfo_file.split("_L001")
    df_allinfo["Sample"] = sample_name
    df_tmp = lookupSubSequencePerFile(df_allinfo)
    df_lookup_long = pd.concat([df_lookup_long, df_tmp])
    print(sample_name, len(df_tmp))
print("ALL", len(df_lookup_long))


# In[ ]:


df_lookup_long.to_excel(outHitsWithMismatch)
print("Wrote", outHitsWithMismatch, "to disk")

