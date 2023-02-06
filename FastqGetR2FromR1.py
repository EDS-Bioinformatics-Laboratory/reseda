import os
import sys
import gzip

# Read fastq files of R1 and store accessions
mydir = "split"
#myfiles = [x for x in os.listdir(mydir) if x.endswith("_R1_001-mapped.fastq.gz")]
myfiles = [x for x in os.listdir(mydir) if x.endswith("_R1_001.fastq.gz")]

# Process in batches of 30
mystart = int(sys.argv[1])
myend = int(sys.argv[2])
myfiles = myfiles[mystart:myend]

accessions = dict()
for f in myfiles:
  # AGTTCTCG_S1_L001_R1_001-mapped.fastq.gz
  print(f)
  barcode = f.split("_")[0]
#  fh = gzip.open(mydir + "/" + f, "rt")
  fh = open(mydir + "/" + f)
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
#      fhOuts[barcode] = gzip.open(mydir + "/" + barcode + "_S1_L001_R2_001-mapped.fastq.gz", "wt")
      fhOuts[barcode] = gzip.open(mydir + "-R2/" + barcode + "_S1_L001_R2_001.fastq.gz", "wt")
    print(descr1, file=fhOuts[barcode])
    print(seq, file=fhOuts[barcode])
    print(descr2, file=fhOuts[barcode])
    print(qual, file=fhOuts[barcode])

for key in fhOuts.keys():
  fhOuts[key].close()
