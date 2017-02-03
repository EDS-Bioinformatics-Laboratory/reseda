from __future__ import print_function
import sys
import os
import json

if __name__ == '__main__':
    jsonFile = "20161216_Rheuma_MiSeqRUN011.json"
    mountdir = "/mnt/immunogenomics/RUNS/run11-20161219-miseq/data/"

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
    index_samples = dict()
    for i in range(len(js["Samples"])):
        sample_name = js["Samples"][i]["Sample_Name"]
        if sample_name in index_samples:               # Check if the sample name is unique
            print("WARNING: multiple samples with same name in "+ run + " " + sample_name)
        index_samples[sample_name] = i

    # Read directory with file paths to the fastq files
    fhs = dict()
    for f in os.listdir(mountdir):
        sample_name_nr = f.split("_L001")
        c = sample_name_nr[0].split("_")
        sample_nr = c.pop(-1)   # last part is the sample number
        sample_name = "_".join(c)  # rest is sample name
        inx = index_samples[sample_name]   # gives an error when sample is not in sample sheet

        # Add the sample number to the json object
        js["Samples"][inx]["Sample_Nr"] = sample_nr

        # Print info to SAMPLES file
        sampleFile = "SAMPLES-"+run+"-"+js["Samples"][inx]["Species"]+"-"+js["Samples"][inx]["Chain"]
        if sampleFile not in fhs:
            fhs[sampleFile] = open(sampleFile, "w")
        print(mountdir+f, file=fhs[sampleFile])

    for sampleFile in fhs:
        fhs[sampleFile].close()
        
    # Write new json file
    jsonNew = jsonFile.replace(".json", "-new.json")
    fhJson = open(jsonNew, "w")
    print(json.dumps(js, indent=4), file=fhJson)
    fhJson.close()