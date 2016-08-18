from __future__ import print_function
import csv
import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def csvToFastq (myfile):
    '''
    Description: read AA.reads file and convert it to fastq
    In: csvfile (str)
    Out: -
    '''
    outfile = myfile + ".fastq.gz"
    fhOut = gzip.open(outfile, "w")
    with open(myfile) as csvfile:
        fh = csv.reader(csvfile)
        for c in fh:
            acc = c[1]
            seq = c[2] + c[3]
            regionmidrun = c[5].split(",")
            run = regionmidrun[-1]
            if run == "142":
                # print(run, acc, seq)
                record = SeqRecord(Seq(seq), id=acc, name=acc)
                record.letter_annotations["phred_quality"] = [40] * len(seq)
                SeqIO.write(record, fhOut, "fastq")
    csvfile.close()
    fhOut.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " AA.reads.csv file(s)")

    for myfile in sys.argv[1:]:
        csvToFastq(myfile)
