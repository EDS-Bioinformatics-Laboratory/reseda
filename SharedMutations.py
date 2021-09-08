import argparse
from Bio import AlignIO

def main(f):

    # Output report files
    outfile = f + ".shared-mutations.txt"
    fhOut = open(outfile, "w")
    mutfile = f + ".mutations.csv"
    fhMut = open(mutfile, "w")

    # Read the MSA made by ClustalO, check which columns belong to the reference sequences and the experimental sequences
    fhIn = open(f)
    aln = AlignIO.read(fhIn, "fasta")
    i = 0
    refs = list()
    exp = list()
    seq_names = list()
    for record in aln:
        #print(record)
        seq_names.append(record.id)
        if record.id.startswith("IG"):
            refs.append(i)
        else:
            exp.append(i)
        i = i + 1
    # print("refs:", refs)
    # print("exp:", exp)
    # print("names:", seq_names)

    # Create a list with nucleotide sequences that are all the same (so we can skip these positions later)
    n = len(aln[:, 0])  # how many sequences?
    same = [n * "c", n * "a", n * "t", n * "g", n * "-"]

    # Create an empty matrix for the mutations
    mutations = dict()
    for e in exp:
        mutations[e] = dict()
        for r in refs:
            mutations[e][r] = list()

    # Go through each position of the MSA and check if nucleotides are different from the germline
    print("Antibody", "Reference", "pos", "ref", "alt", "mut", sep="\t", file=fhMut)
    for i in range(aln.get_alignment_length()):
        column = aln[:, i] # get column

        if column in same: # skip if everything is the same
            continue

        for e in exp: # for every experiment sequence
            for r in refs: # for every reference sequence
                if column[r] != column[e] and column[r] != "-" and column[e] != "-": # check for different nucleotides, skip indels
                    mut = column[r] + str(i+1) + ">" + column[e]
                    mutations[e][r].append(mut)
                    print(seq_names[e], seq_names[r], i+1, column[r], column[e], mut, sep="\t", file=fhMut)

    # Count mutations
    print("Reference", "Antibody", "Mutations", sep="\t", file=fhOut)
    mut_count = dict()
    for r in refs:
        mut_count[r] = dict()
        for e in exp:
            print(seq_names[r], seq_names[e], len(mutations[e][r]), sep="\t", file=fhOut)
            mut_count[r][e] = len(mutations[e][r])

    # Create matrix for common mutations compared to the references
    print("-----", file=fhOut)
    print("Reference", "Antibody 1", "Antibody 2", "Shared mutations", "Perc AB1", "Perc AB2", sep="\t", file=fhOut)
    for r in refs:
        for i in range(len(exp) - 1):
            e1 = exp[i]
            for j in range(i+1, len(exp)):
                e2 = exp[j]
                nr_common_mut = len(set(mutations[e1][r]).intersection(set(mutations[e2][r])))
                perc_e1 = round(100 * nr_common_mut / mut_count[r][e1], 1)
                perc_e2 = round(100 * nr_common_mut / mut_count[r][e2], 1)
                print(seq_names[r], seq_names[e1], seq_names[e2], nr_common_mut, perc_e1, perc_e2, sep="\t", file=fhOut)

    fhOut.close()
    fhMut.close()
    return(mutfile, outfile)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calls mutations from multiple sequence alignment (ClustalO, fasta) and checks for shared mutations')
    parser.add_argument("msa_fasta_files", type=str, nargs='+', help='Path(s) to MSA fasta file(s)')
    args = parser.parse_args()

    for f in args.msa_fasta_files:
        mutfile, outfile = main(f)
        print("Wrote", mutfile, "to disk")
        print("Wrote", outfile, "to disk")
