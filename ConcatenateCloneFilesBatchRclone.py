import json
import sys
import argparse

# Same as ConcatenateCloneFilesBatch.py, but the generated bash script uses rclone
# (like execute-all.sh and rclone-copy-from-webdav.sh) instead of the webdav mount + curl.
#
# Usage:
#   python3 ConcatenateCloneFilesBatchRclone.py -r yyyymmdd-RUNnn-datasheet-new.json -run runNN-yyyymmdd-miseq -o yyyymmdd-reseda > concatenate-clones.sh
#   bash concatenate-clones.sh

# (file pattern after chain_species, prefix of the concatenated file)
FILE_TYPES = [
    ("*.rr.clones_subs.csv", "cdr3-clones-"),
    ("*-all_info.csv.rr.csv", "assign-info-"),
    ("-clones-mut-sites-reassigned.csv", "vjcdr3-clones-mut-"),
]


def executeCmd(cmd):
    '''
    Description: print a command for the generated bash script
    In: command (str)
    Out: -
    '''
    print(cmd)


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
    parser = argparse.ArgumentParser(description='Generates a bash script that downloads (rclone) and concatenates clones files per project and cell type')
    parser.add_argument('-r', '--runinfo', default='yyyymmdd-RUNnn-datasheet-new.json', type=str, help='Sample sheet in json format (default: %(default)s)')
    parser.add_argument('-run', '--run', default='runNN-yyyymmdd-miseq', type=str, help='Run directory on the webdav server (default: %(default)s)')
    parser.add_argument('-o', '--outdir', default='yyyymmdd-reseda', type=str, help='Processing directory of the run (default: %(default)s)')
    parser.add_argument('-p', '--prefix', default='remote', type=str, help='Name of the rclone remote, see: rclone config (default: %(default)s)')
    parser.add_argument('-t', '--tmpdir', default='concatenate-tmp', type=str, help='Local directory for the downloaded files, removed afterwards (default: %(default)s)')
    args = parser.parse_args()

    if args.run == 'runNN-yyyymmdd-miseq' or args.outdir == 'yyyymmdd-reseda' or args.runinfo == 'yyyymmdd-RUNnn-datasheet-new.json':
        parser.print_help()
        sys.exit()

    remotepath = args.prefix + ":RUNS/" + args.run + "/Processing/" + args.outdir + "/Results/final"
    runinfo = args.runinfo
    tmpdir = args.tmpdir

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

    print("#!/bin/bash")
    print("# Remote path: " + remotepath)
    print("#")

    for project, chains_species in projects.items():
        chains_species = list(set(chains_species))
        for chain_specie in chains_species:
            for file_pattern, file_prefix in FILE_TYPES:
                pattern = "*" + chain_specie + file_pattern
                cmd = "rm -rf " + tmpdir + " ; mkdir -p " + tmpdir
                executeCmd(cmd)
                cmd = "rclone copy " + remotepath + " " + tmpdir + " --include '" + pattern + "' --progress --transfers=4 --checkers=8 --retries=3 ; wait"
                executeCmd(cmd)
                cmd = "myfiles=`ls " + tmpdir + "/" + pattern + " 2> /dev/null`"
                executeCmd(cmd)
                cmd = "if [ -n \"$myfiles\" ]; then python3 ConcatenateCloneFiles.py -r " + runinfo + " -n " + project + " -c " + chain_specie + " -pre " + file_prefix + " $myfiles ; else echo 'WARNING: no files found for " + pattern + "' ; fi ; wait"
                executeCmd(cmd)
                cmd = "rm -rf " + tmpdir
                executeCmd(cmd)

                print("#")

    print("echo 'FINISHED'")
