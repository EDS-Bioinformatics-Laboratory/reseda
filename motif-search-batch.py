import sys
import subprocess
import re

#fileList = "hg19-list.txt"
fileList = "AID-list-sep2015.txt"
# fileList = "HLA-list.txt"

try:
    fh = open(fileList)
except:
    sys.exit("cannot open file", fileList)

for myfile in fh:
    myfile = myfile.rstrip()
    if myfile.endswith(".gz"):
        cmd = "gunzip " + myfile
        print cmd
        process = subprocess.Popen(cmd, shell=True, stdout=None, stderr=None, bufsize=1)
        process.wait()
        myfile = re.sub(".gz$", "", myfile)


    myout = open(myfile + ".AID.vcf", 'w')
    # myout = open(myfile + ".primers.txt", 'w')

#    cmd = "python motif-search.py " + myfile
    cmd = "python search-motif-in-genomic-region.py " + myfile
    print cmd
    process = subprocess.Popen(cmd, shell=True, stdout=myout, stderr=None, bufsize=1)
    process.wait()
    
    # Count hits
    # cmd = "cut -d' ' -f 4 " + myfile + ".primers.txt|sort|uniq -c"
    # mycount = open(myfile + ".primers.count.txt", 'w')
    # print cmd
    # process = subprocess.Popen(cmd, shell=True, stdout=mycount, stderr=None, bufsize=1)
    # process.wait()
