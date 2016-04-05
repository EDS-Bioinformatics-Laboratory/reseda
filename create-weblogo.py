from __future__ import print_function
import sys
import os

inFASTA = "weblogo/globin.fasta"
outPDF="globin.pdf"
mycmd = "./weblogo/seqlogo -F PDF -h 20 -w 40 -k 0 -c -f " + inFASTA + " > " + outPDF
print(mycmd)
os.system(mycmd)
print("DONE")
