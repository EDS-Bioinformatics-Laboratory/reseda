from __future__ import print_function
import os
import argparse
import json

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Checks sample sheet (json) and adds sample numbers. Creates a file -new.json')
    parser.add_argument('-r', '--runinfo', default='yyyymmdd-RUNnn-datasheet.json', type=str, help='Sample sheet in json format (default: %(default)s)')
    parser.add_argument('-w', '--webdav', default='/mnt/immunogenomics/RUNS/runNN-yyyymmdd-miseq/data/', type=str, help='Webdav directory (default: %(default)s)')
    args = parser.parse_args()

    if args.webdav == '/mnt/immunogenomics/RUNS/runNN-yyyymmdd-miseq/data/' or args.runinfo == 'yyyymmdd-RUNnn-datasheet.json':
        parser.print_help()
        exit()

    # Read json with parsed sample sheet info (made with MetaData.py)
    try:
        fh = open(args.runinfo)
    except:
        sys.exit("cannot open file: " + args.runinfo)
    text = fh.read()
    js = json.loads(text)
    fh.close()

    # Make a dictionary of the expected samples
    found = dict()
    for sample in js["Samples"]:
        found[sample["Sample_Name"]] = 0

    files = os.listdir(args.webdav)

    check = dict()
    for myfile in files:
        sample, rest = myfile.split("_L001")
        check[sample] = check.get(sample, list())
        check[sample].append(myfile)

    for sample, pairs in check.items():
        sample_name = "_".join(sample.split("_")[:-1])
        if len(pairs) < 2:
            if "R1" in pairs[0]:
                print("grep", sample, "basespace-copy-data.sh | grep R2")
                found[sample_name] = 1
            elif "R2" in pairs[0]:
                print("grep", sample, "basespace-copy-data.sh | grep R1")
                found[sample_name] = 1
        else:
            print(sample, "OK")
            found[sample_name] = 1

    for sample, was_found in found.items():
        if was_found != 1:
            print("grep", sample, "basespace-copy-data.sh")
