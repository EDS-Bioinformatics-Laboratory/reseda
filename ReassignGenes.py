import pandas as pd
import numpy as np
import argparse

def reAssign(df, peptide, threshold):
    '''
    Description: reassign genes with gene with the majority of the reads (or include more genes up to 70% of the reads)
    In: df with all clones, peptide is the CDR3 peptide, threshold is for how many genes need to be included in the description (70% of the reads)
    Out: new V gene name
    '''
    # get all clones for peptide
    df_tmp = df.loc[df['cdr3pep'] == peptide][['cdr3pep','V_sub','acc.nunique']]
    df_tmp = df_tmp.sort_values(by='acc.nunique', ascending=False)

    # calculate total frequency for this peptide
    total_freq = cdr3pep_uniq['acc.nunique'].loc[peptide]

    # store the cumulative sum of the reads
    df_tmp['cumsum'] = df_tmp['acc.nunique'].cumsum()

    # calculate the cumulative fraction of the reads
    df_tmp['cumsum_frac'] = df_tmp['cumsum'] / total_freq

    # select the genes that make up 70% of the reads for this peptide
    include_up_to = [ e > threshold for e in df_tmp['cumsum_frac'].tolist() ].index(True) + 1
    df_select = df_tmp.iloc[0:include_up_to]
    df_select = df_select.sort_values(by='V_sub', ascending=True)

    # concatenate the gene names with a plus sign
    v_genes = list(set(df_select['V_sub'].tolist()))
    v_genes.sort()
    v_gene = "+".join(v_genes)

    # replace v name with new v name
    df.loc[df['cdr3pep'] == peptide, 'V_sub'] = v_gene

    return(df, v_gene)


def mode(array):
    # From: https://stackoverflow.com/questions/10797819/finding-the-mode-of-a-list (user mathwizurd)
    most = max(list(map(array.count, array)))
    return list(set(filter(lambda x: array.count(x) == most, array)))


def min_mode(array):
    return min(pd.Series.mode(array))


def concat_j_gene(genes):
    genes = list(set(genes)) # get unique names
    genes.sort()             # sort alphabetically
    return ",".join(genes)   # concatenate with a comma in between


# def show_mode(mylist):
#     return "|".join(mode(mylist))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reassigns V gene names, creates a clones file')
    parser.add_argument('-c', '--clones', default='SampleName_S1_L001.assembled-MID-IGH_HUMAN-clones-mut-sites.csv', type=str, help='Clone file name, output of MutationAnalysisVJ.py (default: %(default)s)')
    parser.add_argument('-a', '--allinfo', default='final/SampleName_S1_L001.assembled-MID-IGH_HUMAN-all_info.csv', type=str, help='All info file, result from combine-immuno-data.py (default: %(default)s)')
    parser.add_argument('-t', '--threshold', default=0.7, type=float, help='Include V gene names with this cumulative fraction of reads (70 perc, default: %(default)s)')

    args = parser.parse_args()

    if args.allinfo == 'final/SampleName_S1_L001.assembled-MID-IGH_HUMAN-all_info.csv' or args.clones == 'SampleName_S1_L001.assembled-MID-IGH_HUMAN-clones-mut-sites.csv':
        parser.print_help()
        exit()

    # Input
    cloneFile = args.clones
    allinfoFile = args.allinfo
    threshold = args.threshold

    # Output
    outfile = cloneFile.replace("-clones-mut-sites.csv", "-clones-mut-sites-reassigned.csv")
    allinfonew_file = cloneFile.replace("-clones-mut-sites.csv", "-allinfo-filtered.csv")
    reportfile = cloneFile.replace("-clones-mut-sites.csv", "-qual-reassign.log")

    fhOut = open(reportfile, "w")

    # read cloneFile and put it in a dataframe
    df = pd.read_csv(cloneFile, sep="\t", na_values=['None', ''], dtype={'mut.count_x.concat_values': str, 'mut.count_y.concat_values': str, 'nr_sites.concat_values': str})
    print("Nr of entries in clone file", len(df), file=fhOut)
    ###df.head()

    # group by cdr3peptide and count nr of different V genes and reads
    cols = ['cdr3pep', 'V_sub', 'acc.nunique']
    cdr3pep_uniq = df[cols].groupby('cdr3pep').agg({'V_sub': 'nunique', 'acc.nunique': sum})
    print("Nr of unique CDR3s:", len(cdr3pep_uniq), file=fhOut)

    # select CDR3's with more than one V gene assigned
    cdr3pep_uniq = cdr3pep_uniq.loc[cdr3pep_uniq['V_sub'] > 1]
    print("Nr of uniq cdr3s with more than one V gene", len(cdr3pep_uniq), file=fhOut)

    ###cdr3pep_uniq.head()

    ## Loop through all CDR3s and re-assign the V gene ##
    for peptide in cdr3pep_uniq.index:
        (df, new_v_gene) = reAssign(df, peptide, threshold)

    ###df.loc[df['cdr3pep'] == peptide]

    # Group the re-assigned entries
    cols = ['cdr3pep', 'V_sub']
    concat_mode = lambda x: min(mode([float(e) for e in "|".join(x).split("|")]))
    clones = df.groupby(cols).agg({'J_sub': concat_j_gene,'acc.nunique': sum, 'beforeMID.nunique': sum, 'mut.count_x.sum': sum, 'mut.count_x.mean': np.mean, 'mut.count_x.concat_values': concat_mode, 'mut.frac_x.mean': np.mean, 'mut.count_y.sum': sum, 'mut.count_y.mean': np.mean, 'mut.count_y.concat_values': concat_mode, 'mut.frac_y.mean': np.mean, 'nr_sites.sum': sum, 'nr_sites.mean': np.mean, 'nr_sites.concat_values': concat_mode})
    clones = clones.sort_values(by='acc.nunique', ascending=False)
    clones = clones.reset_index()

    ##clones.head()

    ## Check if sum of nr of accessions is the same ##

    print("Sum reads the same?", df['acc.nunique'].sum(), clones['acc.nunique'].sum(), df['acc.nunique'].sum() == clones['acc.nunique'].sum(), file=fhOut)

    print("Sum UMIs the same?", df['beforeMID.nunique'].sum(), clones['beforeMID.nunique'].sum(), df['beforeMID.nunique'].sum() == clones['beforeMID.nunique'].sum(), file=fhOut)

    ## Get nr of unique UMIs from the allinfo file ##

    # Read allinfo file and apply quality filter
    allinfo = pd.read_csv(allinfoFile, sep='\t', na_values=['None', ''])
    print("Reads in allinfo before quality filter:", allinfo['acc'].nunique(), file=fhOut)
    allinfo = allinfo.loc[(allinfo['cdr3_qual_min'] >= 30) & (pd.isna(allinfo['V_sub']) == False) & (pd.isna(allinfo['J_sub']) == False) & ((allinfo['V_flag'] == 0) | (allinfo['V_flag'] == 16)) & ((allinfo['J_flag'] == 0) | (allinfo['J_flag'] == 16))]
    print("Reads in allinfo after quality filter:", allinfo['acc'].nunique(), file=fhOut)

    # Write the allinfo file to disk
    allinfo.to_csv(allinfonew_file, sep="\t", index=False)
    print("Wrote", allinfonew_file, "to disk")

    # Group the original entries by cdr3pep and J-gene
    select = ['cdr3pep', 'V_sub', 'J_sub', 'acc', 'beforeMID', 'cdr3nuc']
    cols = ['cdr3pep']
    clones_orig = allinfo[select].groupby(cols).agg({'beforeMID': 'nunique', 'cdr3nuc': ['nunique', min_mode]})
    clones_orig.columns = ['.'.join(col).strip() for col in clones_orig.columns.values] # Convert multilevel column names to single string column names
    clones_orig = clones_orig.sort_values(by='beforeMID.nunique', ascending=False)

    # Reset index and rename the 'beforeMID' column to 'UMIs'
    clones_orig = clones_orig.reset_index()
    clones_orig = clones_orig.rename(columns={'beforeMID.nunique': 'UMIs'})
    ###clones_orig.head()

    # Merge clones with clones_orig to get the unique number of UMIs
    clones_final = pd.merge(clones, clones_orig, how='inner', left_on=['cdr3pep'], right_on=['cdr3pep'])
    clones_final['UMIs.frac'] = clones_final['UMIs'] / sum(clones_final['UMIs'])
    clones_final = clones_final.sort_values(by='UMIs.frac', ascending=False)
    clones_final = clones_final.rename(columns={"acc.nunique": "freq", "mut.count_x.concat_values": "mut.count_x.mode", "mut.count_y.concat_values": "mut.count_y.mode", "nr_sites.concat_values": "nr_sites.mode"})

    ###clones_final.head()

    print("Nr of clones (after mutation analysis, after V gene assignment, final)", len(clones_orig), len(clones), len(clones_final), file=fhOut)

    # Write the clones to disk
    clones_final.to_csv(outfile, sep='\t', index=False)
    print("Wrote", outfile, "to disk")

    fhOut.close()
