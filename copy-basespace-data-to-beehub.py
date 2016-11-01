from __future__ import print_function
import os

'''
Description: get file names of all fastq files in the basespace directory
In: fill in the sampledir and beehub url below
Out: curl commands on stdout
'''

myurl = "https://beehub.nl/amc-immunogenomics/RUNS/run10-20161031-miseq/data/"

# sampledir = "/data/home/bioinfo/basespace/Runs/RUN06/Samples"
# syscall = os.popen("ls " + sampledir + "/*/Files/Data/Intensities/BaseCalls/*.fastq.gz")
# for line in syscall:
#     line = line.rstrip()
#     print('curl -T "' + line + '" --netrc', myurl, "; wait")

rootdir = "/data/home/bioinfo/basespace/Projects/"
mydirs = ["Paired\ RA", "VDJmouse\ \(2\)", "CEA-CEF", "DNA-RNA", "PsA-SpA"]
for mydir in mydirs:
    syscall = os.popen("ls " + rootdir+mydir+"/Samples/*/Files/*.fastq.gz")

    for line in syscall:
        line = line.rstrip()
        print('curl -T "' + line + '" --netrc', myurl, "; wait")
