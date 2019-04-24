import sys
import argparse
import pandas as pd
import numpy as np
from skbio.diversity import alpha_diversity, get_alpha_diversity_metrics



def calcDiversity(f, sample_col, count_col, diversity_index):
    df = pd.read_csv(f, sep='\t')
    df = df.set_index('cdr3pep')

    # Select columns
    cols = [sample_col, count_col]
    df = df[cols]

    # Make a pivot table
    df_pivot = df.pivot_table(values=count_col, index=sample_col, columns=df.index, aggfunc=np.sum)
    df_pivot = df_pivot.fillna(0)

    # Calculate diversity
    df_diversity = pd.DataFrame({sample_col: df_pivot.index})   # make data frame with sample names
    print(df_diversity.head())
    if diversity_index != 'all':
        df_result = pd.DataFrame({diversity_index: alpha_diversity(diversity_index, df_pivot.values, df_pivot.index)})
        df_result = df_result.reset_index()
        print(df_result.head())
        df_diversity = pd.merge(df_diversity, df_result)
    else:
        for diversity_index in ['berger_parker_d', 'brillouin_d', 'chao1', 'chao1_ci', 'dominance', 'doubles', 'enspie', 'esty_ci', 'fisher_alpha', 'gini_index', 'goods_coverage', 'heip_e', 'kempton_taylor_q', 'margalef', 'mcintosh_d', 'mcintosh_e', 'menhinick', 'observed_otus', 'osd', 'pielou_e', 'robbins', 'shannon', 'simpson', 'simpson_e', 'singles', 'strong']:
            print(diversity_index.upper())
            df_result = pd.DataFrame({diversity_index: alpha_diversity(diversity_index, df_pivot.values, df_pivot.index)})
            df_result = df_result.reset_index()
            print(df_result.head())
            df_diversity = pd.merge(df_diversity, df_result)
    print(df_diversity.head())
    df_diversity.to_csv("diversity-" + f, sep="\t")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate diversity')
    parser.add_argument('-i', '--index', default='simpson', type=str, help='diversity method (default: %(default)s)')
    parser.add_argument('-c', '--count_col', default='UMIs', type=str, help='column with counts (default: %(default)s)')
    parser.add_argument('-s', '--sample_col', default='Sample_Name', type=str, help='column with sample names (default: %(default)s)')
    parser.add_argument("files", type=str, nargs='+', help='Path(s) to clone.csv file(s)')
    args = parser.parse_args()

    # if args.input == "clones.csv":
    #     parser.print_help()
    #     exit()

    for f in args.files:
        calcDiversity(f, args.sample_col, args.count_col, args.index)
