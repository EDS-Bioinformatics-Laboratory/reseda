import pandas as pd
import gzip
from Bio import SeqIO

df = pd.read_csv("bc_celseq2.tsv", sep="\t", header=None)
cellbarcodes = list(df[0])

# Open files to write to for each cell barcode
fh = dict()
for bc in cellbarcodes:
  outfile = "split/" + bc + "_S1_L001_R1_001.fastq.gz"
  fh[bc] = open(outfile, "w")

# Open report file
fh["report"] = open("split/report.txt", "w")

# Open fastq file, look for cell barcode and store sequence in corresponding output file
fastqFile = "AMC-MC-v053_HFCDWGDXL_S1_L001_R1_001.fastq.gz"
for record in SeqIO.parse(gzip.open(fastqFile, "rt"), "fastq"):
  sequence = str(record.seq).upper()
  barcode = sequence[6:14]
  if barcode in cellbarcodes:
    SeqIO.write(record, fh[barcode], "fastq")
    print(record.id, barcode, sep="\t", file=fh["report"])

# Close all the output files
for bc in cellbarcodes:
  fh[bc].close()

fh["report"].close()

