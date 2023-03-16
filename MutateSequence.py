import random
import sys
import gzip
import argparse
from Bio import SeqIO


def changeNucleotide(sequence, pos):
    '''
    Change nucleotide with random other nucleotide and return the change
    '''
    try:
        ref = sequence[pos]
    except:
        sys.exit("ERROR: Something went wrong with position " + pos)
    if ref == "A":
        alt = random.sample(["C", "T", "G"], 1)[0]
    elif ref == "C":
        alt = random.sample(["A", "T", "G"], 1)[0]
    elif ref == "G":
        alt = random.sample(["C", "A", "T"], 1)[0]
    elif ref == "T":
        alt = random.sample(["C", "A", "G"], 1)[0]
    else:
        print("ERROR: unknown nucleotide:", ref)
        exit()
    sequence = sequence[:pos] + alt + sequence[pos + 1:]

    return(sequence, ref, alt)


def mutate_sequence(seq, seq_err):
    '''
    Description: mutates 1 sequence based on given error rate
    In: nucleotide sequence (str), sequence error rate (float)
    Out: mutated sequence (str), list with mutations (position, original nucleotide, new nucleotide)
    '''
    mut_pos = list()
    mutationList = list()
    newseq = seq[:]  # make a copy of the sequence

    for pos in range(len(seq)):
        if random.random() < seq_err:
            mut_pos.append(pos)
    if len(mut_pos) > 0:           # there are mutations
        # mutate positions in sequence
        for pos in mut_pos:
            newseq, ref, alt = changeNucleotide(newseq, pos)
            mutationList.append([str(pos), ref, alt])
    return(newseq, mutationList)


def mutate_all_sequences(f, seq_err):
    '''
    Description: opens fastq file, mutates sequences, stores the result in new fastq file
    In: fastq filename, sequence error rate
    Out: new fastq filename, filename of file with list of mutations
    '''
    fhIn = gzip.open(f, "rt")
    outfile = f.replace(".fastq.gz", "-mutated-" + str(seq_err) + ".fastq.gz")
    mutfile = f.replace(".fastq.gz", "-mutated-" + str(seq_err) + "-list.csv")
    fhOut = gzip.open(outfile, "wt")
    fhMut = open(mutfile, "w")
    print("acc\tpos\tref\talt", file=fhMut)

    for record in SeqIO.parse(fhIn, "fastq"):
        record.seq, mutationList = mutate_sequence(record.seq, seq_err)
        SeqIO.write(record, fhOut, "fastq")
        for mut in mutationList:
            print(record.id, "\t".join(mut), sep="\t", file=fhMut)

    fhIn.close()
    fhOut.close()
    fhMut.close()

    return(outfile, mutfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mutates sequences. Input format: fastq')
    parser.add_argument('-e', '--error', default=0.001, type=float, help='Sequence error rate (default: %(default)s)')
    parser.add_argument("fastq_files", type=str, nargs='+', help='Path(s) to fastq file(s)')
    args = parser.parse_args()

    for f in args.fastq_files:
        print("Mutating", f)
        newfile, mutfile = mutate_all_sequences(f, args.error)
        print("Wrote", newfile, "and", mutfile, "to disk")
