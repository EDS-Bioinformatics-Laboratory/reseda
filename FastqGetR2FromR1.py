import os
import gzip

# Read fastq files of R1 and store accessions
mydir = "align-ALL-V"
myfiles = [x for x in os.listdir(mydir) if x.endswith("_R1_001-mapped.fastq.gz")]

accessions = dict()
for f in myfiles:
  # AGTTCTCG_S1_L001_R1_001-mapped.fastq.gz
  barcode = f.split("_")[0]
  fh = gzip.open(mydir + "/" + f, "rt")
  for descr1 in fh:
    descr1 = descr1.rstrip()
    acc = descr1.split(" ")[0]
    seq = fh.readline()
    descr2 = fh.readline()
    qual = fh.readline()
    accessions[acc] = barcode
  fh.close()

# Read fastq files of R2 and store sequences on disk
fhOuts = dict()
fh = gzip.open("AMC-MC-v053_HFCDWGDXL_S1_L001_R2_001.fastq.gz", "rt")
for descr1 in fh:
  descr1 = descr1.rstrip()
  acc = descr1.split(" ")[0]
  seq = fh.readline()
  seq = seq.rstrip()
  descr2 = fh.readline()
  descr2 = descr2.rstrip()
  qual = fh.readline()
  qual = qual.rstrip()

  if acc in accessions.keys():
    barcode = accessions[acc]
    if barcode not in fhOuts.keys():
      fhOuts[barcode] = gzip.open(mydir + "/" + barcode + "_S1_L001_R2_001-mapped.fastq.gz", "wt")
    print(descr1, file=fhOuts[barcode])
    print(seq, file=fhOuts[barcode])
    print(descr2, file=fhOuts[barcode])
    print(qual, file=fhOuts[barcode])

for key in fhOuts.keys():
  fhOuts[key].close()
