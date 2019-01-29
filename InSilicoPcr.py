from __future__ import print_function
import argparse
import regex
from Bio import SeqIO
import sequences


def buildPattern(primer, mismatches, fwdrev):
    '''
    Description: build a regular expression for the primers
    In: primer (str) (str), mismatches (int), fwdrev (str "fwd" or "rev")
    Out: regular expressions (list of strings)
    '''

    # return empty list when primer is empty
    if primer == "":
        return([])

    primer = primer.upper()
    primers = primer.split(",")

    if mismatches == 0:
        mismatches = ""
    elif mismatches > 0:
        mismatches = "{e<=" + str(mismatches) + "}"
    else:
        raise TypeError('wrong input for mismatches')

    patterns = list()
    for primer in primers:
        if fwdrev == "fwd":
            pattern = "((" + primer + ")" + mismatches + ".+)"
        elif fwdrev == "rev":
            pattern = "(.+" + "(" + primer + ")" + mismatches + ")"
        patterns.append(pattern)

    return(patterns)


def findPrimers(seq_file, seq_format, patterns):
    '''
    Description: find primers (patterns) in all sequences in seq_file
    In: seq_file (str file path), patterns (list of regular expressions)
    Out: print on stdout
    '''
    fh = open(seq_file)
    for record in SeqIO.parse(fh, seq_format):
        seq = str(record.seq).upper()
        comrev_seq = sequences.comrev(seq)
        for pattern in patterns:
            p = regex.compile(pattern, regex.BESTMATCH)
            m = p.search(seq)
            if m is not None:
                (start, end) = m.span(1)
                print(record.id, "+", pattern, start, end, m.group(1))
            m = p.search(comrev_seq)
            if m is not None:
                (start, end) = m.span(1)
                print(record.id, "-", pattern, start, end, m.group(1))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform in-silico PCR')
    parser.add_argument('-fwd', '--forward_primer', default='', type=str, help='Nucleotide sequences separated by comma: CATG,ACGT,ETC (default: %(default)s)')
    # parser.add_argument('-rev', '--reverse_primer', default='', type=str, help='Nucleotide sequence separated by comma: CATG,ACGT,ETC (default: %(default)s)')
    parser.add_argument('-m', '--mismatches', default=0, type=int, help='Allowed nr of mismatches (default: %(default)s)')
    parser.add_argument('-f', '--seq_format', default='fasta', type=str, help='Sequence format, currently only fasta is supported (default: %(default)s)')
    parser.add_argument("seq_files", type=str, nargs='+', help='Path(s) to sequence file(s)')
    args = parser.parse_args()

    patterns = buildPattern(args.forward_primer, args.mismatches, "fwd")

    for seq_file in args.seq_files:
        findPrimers(seq_file, args.seq_format, patterns)
