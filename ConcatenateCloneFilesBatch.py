from __future__ import print_function
# import subprocess
import json
import sys


def executeCmd(cmd):
    '''
    Description: execute a command and check if everything went well
    In: command (str)
    Out: - an exception will be raised when command was not successful
    '''
    print(cmd)
    # cmd = cmd.split()
    # rc = subprocess.call(cmd)
    # if rc != 0:
    #     print("something went wrong with: " + " ".join(cmd))


if __name__ == '__main__':
    mydir = "/mnt/immunogenomics/RUNS/run12-20170127-miseq/results-tbcell/final/correct-mid/"
    runinfo = "/home/barbera/git/tbcell-miseq-pipeline/20170127_Rheuma_MiSeqRUN012.json"

    # Read json file
    try:
        fhJs = open(runinfo, "r")
    except:
        sys.exit("cannot open file")
    text = fhJs.read()
    js = json.loads(text)
    fhJs.close()

    # Make list of all the projects
    projects = list()
    for sample in js["Samples"]:
        projects.append(sample["Sample_Project"])
    projects = list(set(projects))  # make list unique

    for project in projects:
        cmd = "python ConcatenateCloneFiles.py " + runinfo + " " + project + " " + mydir + "*-clones-subs.csv"
        executeCmd(cmd)
        cmd = "python ConcatenateCloneFiles.py " + runinfo + " " + project + " " + mydir + "*.rr.clones_subs.csv"
        executeCmd(cmd)
