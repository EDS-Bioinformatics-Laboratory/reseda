import random
import sys


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
    for pos in range(len(seq)):
        if random.random() < seq_err:
            mut_pos.append(pos)
    if len(mut_pos) > 0:           # there are mutations
        # mutate positions in sequence
        newseq = seq[:] # make a copy of the sequence
        for pos in mut_pos:
            newseq, ref, alt = changeNucleotide(newseq, pos)
            mutationList.append([str(pos), ref, alt])
    print(seq, newseq, mutationList)


if __name__ == '__main__':
    seq_err = 0.1  # 0.001
    seq = "A" * 200
    mutate_sequence(seq, seq_err)
    
