from __future__ import print_function
import sys

''' It doesn't actually roll out the scripts, but it prints the bash code to do so '''

try:
    fh = open("../progress/ip-list")
except:
    sys.exit("cannot open file")

for line in fh:
    line = line.strip()
    # print("scp execute-all.sh bioinfo@" + line + ":~/git/tbcell-miseq-pipeline/")
    # print("scp execute-part.sh bioinfo@" + line + ":~/git/tbcell-miseq-pipeline/")
    # print("scp translate-and-extract-cdr3-1mm.py bioinfo@" + line + ":~/git/tbcell-miseq-pipeline/")
    print("scp setup-and-run.sh bioinfo@" + line + ":~/")

fh.close()
