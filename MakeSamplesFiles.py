from __future__ import print_function
import sys
import os
import json
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Checks sample sheet (json) and adds sample numbers. Creates a file -new.json')
    parser.add_argument('-r', '--runinfo', default='yyyymmdd-RUNnn-datasheet.json', type=str, help='Sample sheet in json format (default: %(default)s)')
    parser.add_argument('-w', '--webdav', default='/mnt/immunogenomics/RUNS/runNN-yyyymmdd-miseq/Data/NameOfDataset_1/Raw/', type=str, help='Webdav directory (default: %(default)s)')
    args = parser.parse_args()

    if args.webdav == '/mnt/immunogenomics/RUNS/runNN-yyyymmdd-miseq/Data/NameOfDataset_1/Raw/' or args.runinfo == 'yyyymmdd-RUNnn-datasheet.json':
        parser.print_help()
        exit()

    jsonFile = args.runinfo   # = "20170709_RUN015_Datasheet_finalV_Selection.json"
    mountdir = args.webdav    # = "/mnt/immunogenomics/RUNS/runNN-yyyymmdd-miseq/Data/NameOfDataset_1/Raw/"

    # Read json with parsed sample sheet info (made with MetaData.py)
    try:
        fh = open(jsonFile)
    except:
        sys.exit("cannot open file: " + jsonFile)
    text = fh.read()
    js = json.loads(text)
    fh.close()

    # run info and sample indices
    run = js["Experiment Name"].lower()
    index_samples = dict()  # To check if there are multiple samples with same name
    file_present = dict()   # To check if files are present for all described samples
    for i in range(len(js["Samples"])):
        sample_name = js["Samples"][i]["Sample_Name"]
        if sample_name in index_samples:               # Check if the sample name is unique
            print("WARNING: multiple samples with same name in " + run + " " + sample_name)
        index_samples[sample_name] = i
        file_present[sample_name] = 0

    # Read directory with file paths to the fastq files
    fhs = dict()
    for f in os.listdir(mountdir):
        if "fastq.gz" not in f:
            continue
        sample_name_nr = f.split("_L001")
        c = sample_name_nr[0].split("_")
        sample_nr = c.pop(-1)   # last part is the sample number
        sample_name = "_".join(c)  # rest is sample name
        try:
            inx = index_samples[sample_name]   # gives an error when sample is not in sample sheet
            file_present[sample_name] = 1
        except:
            print("ERROR:", sample_name, "is not in sample sheet")
            continue

        # Add the sample number to the json object
        js["Samples"][inx]["Sample_Nr"] = sample_nr

        # Check if the project name contains a whitespace. If so, replace it with an underscore
        try:
            if js["Samples"][inx]["Sample_Project"] == "":
                print("WARNING", sample_name, ": 'Sample_Project' is empty")
            js["Samples"][inx]["Sample_Project"] = js["Samples"][inx]["Sample_Project"].replace(" ", "_")
        except:
            print("WARNING", sample_name, ": attribute 'Sample_Project' does not exist")

        # Print info to SAMPLES file
        sampleFile = "SAMPLES-" + run + "-" + js["Samples"][inx]["Species"] + "-" + js["Samples"][inx]["Chain"]
        if sampleFile not in fhs:
            fhs[sampleFile] = open(sampleFile, "w")
        print(mountdir + f, file=fhs[sampleFile])

    # Print all sample names of which there are no files on beehub
    for sample_name, present in file_present.items():
        if present == 0:
            print("Error:", sample_name, "has no files on beehub")

    for sampleFile in fhs:
        fhs[sampleFile].close()

    # Write new json file
    jsonNew = jsonFile.replace(".json", "-new.json")
    fhJson = open(jsonNew, "w")
    print(json.dumps(js, indent=4), file=fhJson)
    fhJson.close()
