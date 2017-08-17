from __future__ import print_function
import sys
from Bio import SeqIO


if __name__ == '__main__':
    imgt_orig = "reference/IGHC_ORIG_human.fasta"
    imgt_selection = "reference/IGHC_CH12_human.fasta"

    try:
        fh = open(imgt_orig)
        fhOut = open(imgt_selection, "w")
    except:
        sys.exit("cannot open files for reading or writing")

    seq_new = dict()
    name_new = dict()
    for record in SeqIO.parse(fh, "fasta"):
        descr_list = record.description.split("|")
        acc = descr_list[0] + "|" + descr_list[1]
        subpart = descr_list[4]
        if "CH1" in subpart or "CH2" in subpart:
            name_new[acc] = name_new.get(acc, acc) + "|" + subpart
            seq_new[acc] = seq_new.get(acc, "") + str(record.seq)

    for acc in name_new:
        print(">" + name_new[acc], file=fhOut)
        print(seq_new[acc], file=fhOut)

    fh.close()
    fhOut.close()
