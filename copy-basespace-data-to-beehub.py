from __future__ import print_function
import os

'''
Description: get file names of all fastq files in the basespace directory
In: fill in the sampledir and beehub url below
Out: curl commands on stdout
'''

myurl = "https://beehub.nl/amc-immunogenomics/RUNS/run12-20170127-miseq/data/"

# sampledir = "/data/home/bioinfo/basespace/Runs/RUN06/Samples"
# syscall = os.popen("ls " + sampledir + "/*/Files/Data/Intensities/BaseCalls/*.fastq.gz")
# for line in syscall:
#     line = line.rstrip()
#     print('curl -T "' + line + '" --netrc', myurl, "; wait")

rootdir = "/data/home/bioinfo/basespace/Projects/"
# mydirs = ["Paired\ RA", "VDJmouse\ \(2\)", "CEA-CEF", "DNA-RNA", "PsA-SpA"]
mydirs = ["VDJmouse\ \(4\)", "Dermatomyositis", "CordBlood\ \(2\)", "CD40L-enrich", "CEA\ \(3\)", "AB-RBF", "AB-RTX", "AB-ETA"]
for mydir in mydirs:
    syscall = os.popen("ls " + rootdir + mydir + "/Samples/*/Files/*.fastq.gz")

    for line in syscall:
        line = line.rstrip()
        print('curl -T "' + line + '" --netrc', myurl, "; wait")
