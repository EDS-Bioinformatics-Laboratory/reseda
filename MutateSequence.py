import random
import sys
import gzip
import argparse
from Bio import SeqIO


def changeNucleotide(sequence, pos):
    """
    Change nucleotide with random other nucleotide and return the change
    """
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
        exit
    sequence = sequence[:pos] + alt + sequence[pos + 1:]

    return(sequence, ref, alt)


def mutate_sequence(seq, seq_err):
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
    print(seq, newseq, mutationList)
    return(newseq)


def mutate_all_sequences(f, seq_err):
    fhIn = gzip.open(f, "rt")
    outfile = f.replace(".fastq.gz", "-mutated-" + str(seq_err) + ".fastq.gz")
    fhOut = gzip.open(outfile, "wt")

    for record in SeqIO.parse(fhIn, "fastq"):
        record.seq = mutate_sequence(record.seq, seq_err)
        SeqIO.write(record, fhOut, "fastq")

    fhIn.close()
    fhOut.close()

    return(outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mutates sequences. Input format: fastq')
    parser.add_argument('-e', '--error', default=0.001, type=float, help='Sequence error rate (default: %(default)s)')
    parser.add_argument("fastq_files", type=str, nargs='+', help='Path(s) to fastq file(s)')
    args = parser.parse_args()

    for f in args.fastq_files:
        print("Mutating", f)
        newfile = mutate_all_sequences(f, args.error)
        print("Wrote", newfile, "to disk")
