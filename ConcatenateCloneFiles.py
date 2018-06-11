from __future__ import print_function
import sys
import json
import re
import argparse


def parseSampleName(myfile):
    '''
    Description: extract the sample name from the file name
    In: file name (str)
    Out: sample name (str), mid (str)
    '''
    path = myfile.split("/")[-1]
    if "_L001.assembled-" in path:
        sample, rest = path.split("_L001.assembled-")
        mid = rest.split("-")[0]
    elif "mutations-per-clone.csv" in path:
        sample_nr_mid = path.replace("-mutations-per-clone.csv", "").split("-")
        sample = "-".join(sample_nr_mid[:-1])
        mid = sample_nr_mid[-1]
    else:
        print("cannot parse sample name and mid from:", path)
        exit()

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
        line = line.replace('"', '')
        line = [sample, mid] + line.strip().split(delimiter)
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
    parser = argparse.ArgumentParser(description='Concatenates clones files for a project and particular cell type')
    parser.add_argument('-r', '--runinfo', default='yyyymmdd-RUNnn-datasheet-new.json', type=str, help='Sample sheet in json format (default: %(default)s)')
    parser.add_argument('-n', '--name_project', default='MYPROJECT', type=str, help='Project name as annotated in the meta data (default: %(default)s)')
    parser.add_argument('-c', '--chain_species', default='IGH_HUMAN', type=str, help='celltype_species (e.g. TRB_HUMAN, IGH_MOUSE) (default: %(default)s)')
    parser.add_argument('-d', '--delimiter', default='\t', type=str, help='Field delimiter (default: %(default)s)')
    parser.add_argument('-pre', '--prefix', default='cdr3-clones-', type=str, help='Prefix of the merged file (default: %(default)s)')
    parser.add_argument("input_files", type=str, nargs='+', help='Path(s) to clones file(s)')
    args = parser.parse_args()

    if args.name_project == 'MYPROJECT' or args.runinfo == 'yyyymmdd-RUNnn-datasheet-new.json':
        parser.print_help()
        exit()

    if len(sys.argv) < 5:
        sys.exit("Usage: concatenate-clone-files.py run-info.json project-name chain-species cdr3-clones- *clones_subs.csv")

    runinfo = args.runinfo
    project_name = args.name_project
    chain_specie = args.chain_species
    prefix = args.prefix
    myfiles = args.input_files
    delimiter = args.delimiter

    print("Project:", project_name)
    print("Chain Specie:", chain_specie)

    # Open file for writing
    try:
        if ".rr." in myfiles[0]:
            fhOut = open(prefix + project_name + "-" + chain_specie + "-after-reassignment.csv", "w")
        else:
            fhOut = open(prefix + project_name + "-" + chain_specie + ".csv", "w")
    except:
        sys.exit("cannot write to file")

    # Read header of the first file and write to disk
    try:
        fh = open(myfiles[0])
    except:
        sys.exit("cannot open file " + myfiles[0])

    header = fh.readline()
    header = header.replace('"', '')
    header = ["Sample", "MID"] + header.strip().split(delimiter)
    print("\t".join(header), file=fhOut)
    fh.close()

    # Determine which samples to include
    samples_of_project = getSamplesOfProject(project_name, runinfo)

    # Read content of all files and write to fhOut if sample belongs to project
    p = re.compile("_S\d+$")   # Ends with _S33
    for myfile in myfiles:
        sample_name, mid = parseSampleName(myfile)
        sample_name = p.sub("", sample_name)
        if sample_name in samples_of_project:
            printContent(myfile, fhOut)

    fhOut.close()

    print("DONE")
