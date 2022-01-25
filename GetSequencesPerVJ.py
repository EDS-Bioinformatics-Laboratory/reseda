import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def readAllInfo(f):
    samplename, rest = f.split("/")[-1].split("_L001")
    df = pd.read_csv(f, sep="\t")
    cols = ['V_gene', 'J_gene', 'cdr3nuc', 'acc', 'seq']
    df = df[cols]
    df["Sample"] = samplename
    return(df)

def selectVJ(df, option_vgene, option_jgene):
    # Which VJ combination occurs most frequent?
    df_vj_count = df.groupby(['V_gene', 'J_gene']).agg({'Sample': 'nunique', 'acc': 'nunique'})
    df_vj_count = df_vj_count.sort_values(by=['Sample','acc'], ascending=False)
    df_vj_count = df_vj_count.reset_index()
    df_vj_count.head()

    # Select entries with specified V and J gene (incl allele) or select the top VJ combination
    if option_vgene != "top" and option_jgene != "top": # Were both the V and J genes given as an option? Select those
        v_top = option_vgene
        j_top = option_jgene
    else:                                      # Else, select the VJ combination that occurs often
        v_top = df_vj_count['V_gene'][0]
        j_top = df_vj_count['J_gene'][0]

    # Retrieve sequences for VJ combination and make these sequences unique for further analysis
    concatenate = lambda x: "|".join(list(set(x)))
    df_selection = df[(df["V_gene"] == v_top) & (df["J_gene"] == j_top)]
    df_selection = df_selection.groupby("seq").agg({'acc': concatenate, 'Sample': concatenate, 'cdr3nuc': concatenate}).reset_index()
    print("Selected entries:", len(df_selection))

    return(df_selection, v_top, j_top)

def writeFasta(df_selection, v_top, j_top, option_outdir, option_project_name, option_comrev):
    # Open files for writing the fasta sequences
    out_fastaV = option_outdir + option_project_name + ".V.fasta"
    out_fastaJ = option_outdir + option_project_name + ".J.fasta"
    fhOutV = open(out_fastaV, "w")
    fhOutJ = open(out_fastaJ, "w")

    # Retrieve the reference V and J sequences
    v_seq, j_seq = "", ""
    for record in SeqIO.parse(open("./reference/IGHV_human.fasta"), "fasta"):
        if v_top in record.id:
            v_seq = str(record.seq).upper()
    for record in SeqIO.parse(open("./reference/IGHJ_human.fasta"), "fasta"):
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

    # Function reverse complement
    comrev = lambda s: str(Seq(s).reverse_complement()).upper()

    # Write the sequences in fasta format to disk
    for i in df_selection.index:
        print(">" + df_selection.iloc[i]['acc'] + "|" + df_selection.iloc[i]['Sample'] + "|" + df_selection.iloc[i]['cdr3nuc'], file=fhOutV)
        if option_comrev == 1:
            print(comrev(df_selection.iloc[i]['seq']), file=fhOutV)
        else:
            print(df_selection.iloc[i]['seq'], file=fhOutV)

        print(">" + df_selection.iloc[i]['acc'] + "|" + df_selection.iloc[i]['Sample'] + "|" + df_selection.iloc[i]['cdr3nuc'], file=fhOutJ)
        if option_comrev == 1:
            print(comrev(df_selection.iloc[i]['seq']), file=fhOutJ)
        else:
            print(df_selection.iloc[i]['seq'], file=fhOutJ)

    fhOutV.close()
    fhOutJ.close()

    return(out_fastaV, out_fastaJ)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reads allinfo table, select entries for the top VJ or specified VJ combination, writes sequences to disk')
    parser.add_argument('-v', '--vgene', default='top', type=str, help='V gene, including the allele, e.g. IGHV3-23*01 (default: %(default)s)')
    parser.add_argument('-j', '--jgene', default='top', type=str, help='J gene, including the allele, e.g. IGHJ4*02 (default: %(default)s)')
    parser.add_argument('-cr', '--comrev', default=0, type=int, help='Make experiment reverse_complement (0: no, 1: yes) (default: %(default)s)')
    parser.add_argument('-p', '--project', default='myproject', type=str, help='Project name, this is the prefix of the output fasta files (default: %(default)s)')
    parser.add_argument('-o', '--outdir', default='./', type=str, help='Output directory for fasta files, with trailing / (default: %(default)s)')
    parser.add_argument("allinfo_files", type=str, nargs='+', help="Path(s) to allinfo file(s). Tables should contain the columns 'V_gene', 'J_gene', 'cdr3nuc', 'acc' and 'seq')")

    args = parser.parse_args()

    if args.project == 'myproject':
        print("Please specify a project name / prefix for your output files")
        parser.print_help()
        exit()

    # Read the all_info files (tables should contain the columns 'V_gene', 'J_gene', 'cdr3nuc', 'acc' and 'seq')
    df = readAllInfo(args.allinfo_files[0])
    for f in args.allinfo_files[1:]:
        df = pd.concat([df, readAllInfo(f)])

    # Make a selection with the specified V and J genes or select the top VJ combination
    df_selection, v_top, j_top = selectVJ(df, args.vgene, args.jgene)

    # Write the sequences to disk in fasta format
    out_fastaV, out_fastaJ = writeFasta(df_selection, v_top, j_top, args.outdir, args.project, args.comrev)
    print("Wrote", out_fastaV, "to disk")
    print("Wrote", out_fastaJ, "to disk")
