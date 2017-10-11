from __future__ import print_function
import sys
import os
import gzip
from Bio import SeqIO
import argparse


def splitFastq(f, seq_length):
    '''
    Description: split fastq file based on sequence length
    In: f (file name, str), seq_length (int)
    Out: new files will be written to disk, count_short (int), count_long (int), count_total (int)
    '''
    try:
        fh = gzip.open(f, "rb")
    except:
        sys.exit("cannot open file:" + f)

    fastq_short = f.replace(".fastq.gz", ".short." + str(seq_length) + ".fastq.gz")
    fastq_long = f.replace(".fastq.gz", ".long." + str(seq_length) + ".fastq.gz")
    try:
        fh_short = gzip.open(fastq_short, "w")
        fh_long = gzip.open(fastq_long, "w")
    except:
        sys.exit("cannot write to file: " + fastq_short + " or " + fastq_long)

    count_total = 0
    count_short = 0
    count_long = 0
    for record in SeqIO.parse(fh, "fastq"):
        if len(record.seq) < seq_length:
            SeqIO.write(record, fh_short, "fastq")
            count_short += 1
        else:
            SeqIO.write(record, fh_long, "fastq")
            count_long += 1
        count_total += 1

    fh.close()
    fh_short.close()
    fh_long.close()
    print("Wrote", fastq_short, "to disk")
    print("Wrote", fastq_long, "to disk")
    return(count_short, count_long, count_total)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Split fastq files based on sequence length. Threshold will be used as minimum for the longer sequences (greater than or equal to)')
    parser.add_argument('-l', '--length', default=270, type=int, help='Threshold: sequence length (default: %(default)s)')
    parser.add_argument("fastq_files", type=str, nargs='+', help='Path(s) to fastq file(s)')
    args = parser.parse_args()

    fh_report = open("report-SEQLENGTH.txt", "w")
    print("File Total Short Long Perc_short Perc_long Threshold", file=fh_report)
    for fastq_file in args.fastq_files:
        (count_short, count_long, count_total) = splitFastq(fastq_file, args.length)
        perc_short = 100 * count_short / float(count_total)
        perc_long = 100 * count_long / float(count_total)
        print(fastq_file, count_total, count_short, count_long, perc_short, perc_long, args.length, file=fh_report)
    fh_report.close()
    print("Wrote report-SEQLENGTH.txt to disk")
