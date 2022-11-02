import sys

def renameFile(f):
  sample, pair = f.split("_")
  pair = pair.replace(".fastq.gz", "")
  new_f = sample + "_S1_L001_R" + pair + "_001.fastq.gz"
  print("mv", f, new_f)

for f in sys.argv[1:]:
  renameFile(f)

