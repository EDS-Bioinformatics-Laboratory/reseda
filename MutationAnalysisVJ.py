import pandas as pd
import numpy as np
import argparse


def readData(allinfo_file, v_file, j_file):
    # ## Read files
    allinfo = pd.read_csv(allinfo_file, sep='\t', na_values=['None', ''])

    # Replace 'None' with 0 for the nr_sites column
    allinfo['nr_sites'] = allinfo['nr_sites'].fillna(0)

    v = pd.read_csv(v_file, sep=' ', na_values=['None', ''])
    j = pd.read_csv(j_file, sep=' ', na_values=['None', ''])

    # clean up the gene names
    clean_name = lambda x: x.split("|")[1]
    v['gene'] = [g for g in map(clean_name, v['gene'])]
    j['gene'] = [g for g in map(clean_name, j['gene'])]

    return(allinfo, v, j)

def combineDataFrames(allinfo, v, j):
    # ## Combine files
    df = pd.merge(allinfo, v, how='left', left_on=['acc','V_gene'], right_on=['acc','gene'])
    df = pd.merge(df, j, how='left', left_on=['acc','J_gene'], right_on=['acc','gene'])
    return(df)


def filterData(df):
    # ## Filter data
    df = df.loc[(df['cdr3_qual_min'] >= 30) & (pd.isna(df['V_sub']) == False) & (pd.isna(df['J_sub']) == False) & ((df['V_flag'] == 0) | (df['V_flag'] == 16)) & ((df['J_flag'] == 0) | (df['J_flag'] == 16))]

    # Remove entries where the V and J alignments overlap each other
    df = df.drop(df.loc[(df['start.pos_y']>df['start.pos_x']) & (df['start.pos_y']<df['end.pos_x']) | (df['end.pos_y']>df['start.pos_x']) & (df['end.pos_y']<df['end.pos_x'])].index)

    # Select the alignment with the longest alignment length (V gene)
    longest_alignment = df.groupby('acc').agg({'align.length_x': max})
    longest_alignment = longest_alignment.reset_index()

    df = pd.merge(df, longest_alignment, how='inner', left_on=['acc','align.length_x'], right_on=['acc','align.length_x'])

    return(df)


# Concatenate mutation count values with a '|' between the values
def concat_values(x):
    return "|".join([str(e) for e in x])


def allinfoToClones(df,allinfo_file):
    # ## Group data per clone (CDR3pep)

    # ## Output
    outfile = allinfo_file.replace("-all_info.csv", "-clones-mut-sites.csv")

    clones = df.groupby(['cdr3pep','V_sub','J_sub']).agg({'acc': 'nunique', 'beforeMID': 'nunique', 'mut.count_x': [sum, np.mean, concat_values], 'mut.frac_x': np.mean, 'mut.count_y': [sum, np.mean, concat_values], 'mut.frac_y': np.mean, 'nr_sites': [sum, np.mean, np.median, concat_values]})
    clones = clones.sort_values(by=('acc','nunique'), ascending=False)

    # Convert multilevel column names to single string column names
    clones.columns = ['.'.join(col).strip() for col in clones.columns.values]

    clones.to_csv(outfile, sep='\t')
    print("Wrote", outfile, "to disk")

    return(clones)


def main(allinfo_file, v_file, j_file):
    (allinfo, v, j) = readData(allinfo_file, v_file, j_file)
    df = combineDataFrames(allinfo, v, j)
    df.to_csv("try.csv", sep='\t')

    df = filterData(df)
    outAllinfoFiltered = allinfo_file.replace("-all_info.csv", "-allinfo-filtered-mut.csv")
    df.to_csv(outAllinfoFiltered)
    print("Wrote", outAllinfoFiltered, "to disk")

    clones = allinfoToClones(df, allinfo_file)
    print("FINISHED")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine allinfo with mutations in V and J, creates a clones file')
    parser.add_argument('-a', '--allinfo', default='final/SampleName_S1_L001.assembled-MID-IGH_HUMAN-all_info.csv', type=str, help='All info file, result from combine-immuno-data.py (default: %(default)s)')
    parser.add_argument('-v', '--v', default='SampleName_S1_L001.assembled-MID-IGHV_human-e-clean.sam.mut.txt', type=str, help='Cleaned sam file for V gene assignment (default: %(default)s)')
    parser.add_argument('-j', '--j', default='SampleName_S1_L001.assembled-MID-IGHJ_human-e-clean.sam.mut.txt', type=str, help='Cleaned sam file for J gene assignment (default: %(default)s)')

    args = parser.parse_args()

    if args.allinfo == 'final/SampleName_S1_L001.assembled-MID-IGH_HUMAN-all_info.csv' or args.v == 'SampleName_S1_L001.assembled-MID-IGHV_human-e-clean.sam.mut.txt' or args.j == 'SampleName_S1_L001.assembled-MID-IGHJ_human-e-clean.sam.mut.txt':
        parser.print_help()
        exit()

    main(args.allinfo, args.v, args.j)
