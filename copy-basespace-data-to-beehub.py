from __future__ import print_function
import os

'''
Description: get file names of all fastq files in the basespace directory
In: fill in the sampledir and beehub url below
Out: curl commands on stdout
'''

myurl = "https://beehub.nl/amc-immunogenomics/RUNS/run16-20170808-miseq/data/"

# sampledir = "/data/home/bioinfo/basespace/Runs/RUN06/Samples"
# syscall = os.popen("ls " + sampledir + "/*/Files/Data/Intensities/BaseCalls/*.fastq.gz")
# for line in syscall:
#     line = line.rstrip()
#     print('curl -T "' + line + '" --netrc', myurl, "; wait")

rootdir = "/data/home/bioinfo/basespace/Projects/"
# mydirs = ["Paired\ RA", "VDJmouse\ \(2\)", "CEA-CEF", "DNA-RNA", "PsA-SpA"]
# mydirs = ["VDJmouse\ \(4\)", "Dermatomyositis", "CordBlood\ \(2\)", "CD40L-enrich", "CEA\ \(3\)", "AB-RBF", "AB-RTX", "AB-ETA"]
# mydirs = ["AB-RBF-35707717", "CEA-CEF-35710715"]
mydirs = ["ABI-RTX\ \(2\)", "Cordblood\ \(2\)", "isolationtest", "CEA-CEF\ \(4\)", "ABI-MS\ \(2\)", "SolMem", "TEC", "C1", "CEA-Ada\ \(2\)", "Vasculitis\ \(2\)"]
for mydir in mydirs:
    syscall = os.popen("ls " + rootdir + mydir + "/Samples/*/Files/*.fastq.gz")
    # syscall = os.popen("ls " + rootdir + mydir + "/*/*.fastq.gz")

    for line in syscall:
        line = line.rstrip()
        print('curl -T "' + line + '" --netrc', myurl, "; wait")
