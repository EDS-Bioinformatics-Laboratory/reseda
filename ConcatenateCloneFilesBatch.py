from __future__ import print_function
# import subprocess
import json
import sys
import argparse

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
    parser = argparse.ArgumentParser(description='Concatenates clones files per project and cell type')
    parser.add_argument('-r', '--runinfo', default='yyyymmdd-RUNnn-datasheet-new.json', type=str, help='Sample sheet in json format (default: %(default)s)')
    parser.add_argument('-w', '--webdav', default='/mnt/immunogenomics/RUNS/runNN-yyyymmdd-miseq/results-tbcell/final/', type=str, help='Webdav directory (default: %(default)s)')
    args = parser.parse_args()

    if args.webdav == '/mnt/immunogenomics/RUNS/runNN-yyyymmdd-miseq/results-tbcell/final/' or args.runinfo == 'yyyymmdd-RUNnn-datasheet-new.json':
        parser.print_help()
        exit()

    mydir = args.webdav
    runinfo = args.runinfo

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
            cmd = "python ConcatenateCloneFiles.py -r " + runinfo + " -n " + project + " -c " + chain_specie + " -pre " + "cdr3-clones-" + " " + mydir + "*" + chain_specie + "*.rr.clones_subs.csv"
            executeCmd(cmd)
            cmd = "python ConcatenateCloneFiles.py -r " + runinfo + " -n " + project + " -c " + chain_specie + " -pre " + "assign-info-" + " " + mydir + "*" + chain_specie + "*-all_info.csv.rr.csv"
            executeCmd(cmd)
            cmd = "python ConcatenateCloneFiles.py -r " + runinfo + " -n " + project + " -c " + chain_specie + " -pre " + "vjcdr3-clones-mut-" + " " + mydir + "*-clones-mut-sites-reassigned.csv"
            executeCmd(cmd)

    print("echo 'FINISHED'")
