from __future__ import print_function
import os

'''
Description: get file names of all fastq files in the basespace directory
In: fill in the sampledir and beehub url below
Out: curl commands on stdout
'''

sampledir = "/data/home/bioinfo/basespace/Runs/RUN06/Samples"
myurl = "https://beehub.nl/amc-immunogenomics/RUNS/run06-20160306-miseq/data/"

syscall = os.popen("ls " + sampledir + "/*/Files/Data/Intensities/BaseCalls/*.fastq.gz")

for line in syscall:
    line = line.rstrip()
    print('curl -T "' + line + '" --netrc', myurl, "; wait")
