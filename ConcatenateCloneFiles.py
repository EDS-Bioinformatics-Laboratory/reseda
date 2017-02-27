from __future__ import print_function
import sys
import json
import re


def parseSampleName(myfile):
    '''
    Description: extract the sample name from the file name
    In: file name (str)
    Out: sample name (str)
    '''
    path = myfile.split("/")[-1]
    sample, rest = path.split("_L001.assembled-")
    mid = rest.split("-")[0]
    return(sample, mid)


def printContent(myfile, fhOut):
    '''
    Description: open file, skip header and print to file
    '''
    try:
        fh = open(myfile)
    except:
        sys.exit("cannot open file " + myfile)

    sample, mid = parseSampleName(myfile)

    fh.readline()  # skip header
    for line in fh:
        line = [sample, mid] + line.strip().split("\t")
        print("\t".join(line), file=fhOut)

    fh.close()


def getSamplesOfProject(project_name, runinfo):
    '''
    Description: Go through the run information and get the samples of project
    In: project_name (str), runinfo (path to runXX.json)
    Out: samples (list of sample names)
    '''

    # Read json file
    try:
        fhJs = open(runinfo, "r")
    except:
        sys.exit("cannot open file")
    text = fhJs.read()
    js = json.loads(text)
    fhJs.close()

    samples = list()
    for sample in js["Samples"]:
        if sample["Sample_Project"] == project_name:
            samples.append(sample["Sample_Name"])

    return(samples)


if __name__ == '__main__':

    if len(sys.argv) < 2:
        sys.exit("Usage: concatenate-clone-files.py run-info.json project-name *clones_subs.csv")

    runinfo = sys.argv[1]
    project_name = sys.argv[2]
    myfiles = sys.argv[3:]

    # Open file for writing
    try:
        if ".rr." in myfiles[0]:
            fhOut = open("run-clones_subs-" + project_name + "-after-reassignment.csv", "w")
        else:
            fhOut = open("run-clones_subs-" + project_name + ".csv", "w")
    except:
        sys.exit("cannot write to file")

    # Read header of the first file and write to disk
    try:
        fh = open(myfiles[0])
    except:
        sys.exit("cannot open file " + myfiles[0])

    header = fh.readline()
    header = ["Sample", "MID"] + header.strip().split("\t")
    print("\t".join(header), file=fhOut)
    fh.close()

    # Determine which samples to include
    samples_of_project = getSamplesOfProject(project_name, runinfo)

    # Read content of all files and write to fhOut if sample belongs to project
    p = re.compile("_S\d+$")   # Ends with _S33
    for myfile in myfiles:
        path = myfile.split("/")[-1]
        sample_name, rest = path.split("_L001.assembled-")
        sample_name = p.sub("", sample_name)
        if sample_name in samples_of_project:
            printContent(myfile, fhOut)

    fhOut.close()

    print("DONE")
