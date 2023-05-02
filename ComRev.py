import sys
from Bio import SeqIO
from Bio import Seq

for f in sys.argv[1:]:
    fhOut = open(f + "-comrev.fa", "w")

    for record in SeqIO.parse(f, "fasta"):
        rec_comrev = record.reverse_complement()
        rec_comrev.id = record.id
        rec_comrev.description = record.description
        SeqIO.write(rec_comrev, fhOut, "fasta")
    fhOut.close()

