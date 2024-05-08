import os
import sys
import json
import subprocess
from CleanupDirectory import cleanUp

def readJsonFile(f):
    # Read json file
    try:
        fh = open(f)
    except:
        sys.exit("cannot open file: " + f)
    text = fh.read()
    js = json.loads(text)
    fh.close()
    return(js)

def runJob(js):
    '''
    Description: Prepare and run job
    In: js (json)
    Out: return-code (int)
    '''

    # Write SAMPLES file
    try:
        fhOut = open("SAMPLES", "w")
    except:
        sys.exit("cannot write SAMPLES file")
    for sample in js["samples"]:
        print(sample, file=fhOut)
    fhOut.close()

    # Start analysis
    cmd = ["./execute-all.sh", "-r", js["run"], "-m", js["mids"], "-org", js["organism"], "-cell", js["cell"],
         "-celltype", js["celltype"], "-mm", str(js["mismatches"]), "-s", str(js["seqlength"]), "-cregion",
         js["cregion"], "-p", js["protocol"], "-o", js["outdir"], "-b", js["barcodes"], "-u", js["umis"]]
    print(cmd)
    rc = subprocess.call(cmd)
    return(rc)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage:", sys.argv[0], "tokens_directory/")
        exit()

    tokens_dir = sys.argv[1]
    json_files = [x for x in os.listdir(tokens_dir) if x.endswith(".json")]

    # Read json file with all the parameters, execute RESEDA, and clean-up
    for f in json_files:
        js = readJsonFile(tokens_dir + f) # get parameters
        runJob(js)           # execute RESEDA
        cleanUp()            # remove output files

