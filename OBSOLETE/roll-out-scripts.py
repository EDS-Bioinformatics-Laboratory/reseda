from __future__ import print_function
import sys

''' It doesn't actually roll out the scripts, but it prints the bash code to do so '''

try:
    fh = open("../progress/ip-list")
except:
    sys.exit("cannot open file")

for line in fh:
    line = line.strip()
    print("scp setup-and-run.sh bioinfo@" + line + ":~/")

fh.close()
