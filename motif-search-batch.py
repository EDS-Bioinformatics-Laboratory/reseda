import sys
import subprocess
import re

if len(sys.argv)<2:
    sys.exit("Usage: motif-search-batch.py fastq-file(s)")

fileList = sys.argv[1:]

for myfile in fileList:
    myfile = myfile.rstrip()
    # if myfile.endswith(".gz"):
    #     cmd = "gunzip " + myfile
    #     print cmd
    #     process = subprocess.Popen(cmd, shell=True, stdout=None, stderr=None, bufsize=1)
    #     process.wait()
    #     myfile = re.sub(".gz$", "", myfile)


#    myout = open(myfile + ".AID.vcf", 'w')
    myout = open(myfile + ".primers.txt", 'w')

    cmd = "python motif-search.py " + myfile
    # cmd = "python search-motif-in-genomic-region.py " + myfile
    print cmd
    process = subprocess.Popen(cmd, shell=True, stdout=myout, stderr=None, bufsize=1)
    process.wait()
    
    # Count hits
    cmd = "cut -d' ' -f 4 " + myfile + ".primers.txt|sort|uniq -c|sort -nr"
    mycount = open(myfile + ".primers.count.txt", 'w')
    print cmd
    process = subprocess.Popen(cmd, shell=True, stdout=mycount, stderr=None, bufsize=1)
    process.wait()
