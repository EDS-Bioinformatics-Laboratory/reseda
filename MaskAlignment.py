from __future__ import print_function
import sys
import argparse
import gzip
import regex
from Bio import SeqIO
import MutationsFromSam


def readFastq(f):
    try:
        fh = gzip.open(f, "r")
    except:
        sys.exit("cannot open file:" + f)

    fastq = dict()
    for record in SeqIO.parse(fh, "fastq"):
        acc = record.id
        sequence = str(record.seq).upper()
        fastq[acc] = sequence
    fh.close()
    return(fastq)


def maskAlignment(fastq, sam, letter):
    try:
        fh = open(sam)
    except:
        sys.exit("cannot open file: " + sam)

    for line in fh:
        # Skip header
        if line.startswith("@"):
            continue

        # Read alignment
        line = line.rstrip()
        line = line.split()
        try:
            cigar = MutationsFromSam.parseCigar(line[5])
        except:
            continue
        acc = line[0]
        seq = line[9]
        if len(cigar) > 1:
            print("WARNING: multiple parts aligned:", acc, line[5], cigar, "\n", seq, "\n-----------------------------------------------------------------------")
        for (start, end) in cigar:
            fastq[acc] = fastq[acc][:start] + letter * (end - start) + fastq[acc][end:]
    return(fastq)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Masks aligned areas from sequences')
    parser.add_argument('fastq', type=str, help='fastq.gz file')
    parser.add_argument('v_sam', type=str, help='Alignment file V gene')
    parser.add_argument('j_sam', type=str, help='Alignment file J gene')
    parser.add_argument('-c', '--c_sam', type=str, help='Alignment file C gene (optional)')
    args = parser.parse_args()

    # Read fastq and put it in a dictionary
    fastq = readFastq(args.fastq)

    # Read V alignment file and mask sequences based on the alignment (the cigar string)
    fastq = maskAlignment(fastq, args.v_sam, "v")

    # Read J alignment file and mask sequences based on the alignment (the cigar string)
    fastq = maskAlignment(fastq, args.j_sam, "j")

    # Read C alignment file and mask sequences based on the alignment (the cigar string)
    if args.c_sam is not None:
        fastq = maskAlignment(fastq, args.c_sam, "c")

    # Find part between v and j, print everything
    print("acc orientation np_additions seq")
    p_plus = regex.compile("v(.*)j")
    p_min = regex.compile("j(.*)v")
    for acc, seq in fastq.items():
        np_additions = ""
        direction = "?"
        m_plus = p_plus.search(seq)
        m_min = p_min.search(seq)
        if m_plus is not None:
            np_additions = m_plus.group(1)
            direction = "+"
        elif m_min is not None:
            np_additions = m_min.group(1)
            direction = "-"
        np_additions = np_additions.replace("v", "")
        np_additions = np_additions.replace("j", "")
        print(acc, direction, np_additions, seq)
