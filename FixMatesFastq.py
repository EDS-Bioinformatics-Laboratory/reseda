# This script assumes that the read1 and read2 fastq files are sorted on name (SortFastq.py)
# Usage: python FixMatesFastq.py read1.fastq.gz read2.fastq.gz

import sys
import gzip

def read_fastq(f):
  header = ''
  data = ''

  allData = dict()

  fh = gzip.open(f, "rt")
  for line in fh:
      if line[0] == "@":
          if data != '':
              allData[header] = data

          header = line.strip()[1:].split(" ")[0]
          data = line
      else:
          data += line

  allData[header] = data
  fh.close()
  return(allData)

fastq1 = sys.argv[1]
fastq2 = sys.argv[2]

data1 = read_fastq(fastq1)
data2 = read_fastq(fastq2)

fqOut1 = fastq1.replace(".fastq.gz", ".fixmates.fastq.gz")
fqOut2 = fastq2.replace(".fastq.gz", ".fixmates.fastq.gz")

accs1 = set(data1.keys())
accs2 = set(data2.keys())
intersect = accs1.intersection(accs2)

fhOut1 = gzip.open(fqOut1, "wt")
fhOut2 = gzip.open(fqOut2, "wt")

for acc in intersect:
  print(data1[acc], end='', file=fhOut1)
  print(data2[acc], end='', file=fhOut2)

fhOut1.close()
fhOut2.close()

