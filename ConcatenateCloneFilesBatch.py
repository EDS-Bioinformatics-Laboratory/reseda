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


def lookupChain(chain):
    if chain.lower() == "tcrb":
        chain = "TRB"
    elif chain.lower() == "tcra":
        chain = "TRA"
    elif chain.lower() == "bcrh" or chain.lower() == "igh" or chain.lower() == "bcr-ig":
        chain = "IGH"
    elif chain.lower() == "bcrl" or chain.lower() == "igl":
        chain = "IGL"
    elif chain.lower() == "bcrk" or chain.lower() == "igk":
        chain = "IGK"
    return(chain)


if __name__ == '__main__':
    mydir = "/mnt/immunogenomics/RUNS/run17-20170821-miseq/results-tbcell/final/correct-mid/"
    runinfo = "20170821-RUN17-datasheet-new.json"

    # Read json file
    try:
        fhJs = open(runinfo, "r")
    except:
        sys.exit("cannot open file")
    text = fhJs.read()
    js = json.loads(text)
    fhJs.close()

    # Make list of all the projects
    projects = dict()
    for sample in js["Samples"]:
        projects[sample["Sample_Project"]] = projects.get(sample["Sample_Project"], list())
        species = sample["Species"].upper()
        chain = lookupChain(sample["Chain"])
        projects[sample["Sample_Project"]].append(chain + "_" + species)

    for project, chains_species in projects.items():
        chains_species = list(set(chains_species))
        for chain_specie in chains_species:
            #cmd = "python ConcatenateCloneFiles.py " + runinfo + " " + project + " " + chain_specie + " " + mydir + "*" + chain_specie + "*-clones-subs.csv"
            #executeCmd(cmd)
            cmd = "python ConcatenateCloneFiles.py " + runinfo + " " + project + " " + chain_specie + " " + mydir + "*" + chain_specie + "*.rr.clones_subs.csv"
            executeCmd(cmd)
