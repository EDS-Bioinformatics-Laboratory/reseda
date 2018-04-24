from __future__ import print_function
import os
import argparse
import json
import subprocess

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Checks sample sheet (json) and adds sample numbers. Creates a file -new.json')
    parser.add_argument('-i', '--info', default='yyyymmdd-RUNnn-datasheet.json', type=str, help='Sample sheet in json format (default: %(default)s)')
    parser.add_argument('-r', '--run', default='runNN-yyyymmdd-miseq', type=str, help='Run name (default: %(default)s)')
    parser.add_argument('-wm', '--webdav_mount', default='/mnt/immunogenomics/RUNS/', type=str, help='Webdav mounted directory (default: %(default)s)')
    parser.add_argument('-wu', '--webdav_url', default='https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/', type=str, help='Webdav mounted directory (default: %(default)s)')

    args = parser.parse_args()

    if args.run == 'runNN-yyyymmdd-miseq' or args.info == 'yyyymmdd-RUNnn-datasheet.json':
        parser.print_help()
        exit()

    webdav = args.webdav_mount + args.run + "/data/"

    # Read json with parsed sample sheet info (made with MetaData.py)
    try:
        fh = open(args.info)
    except:
        sys.exit("cannot open file: " + args.info)
    text = fh.read()
    js = json.loads(text)
    fh.close()

    # Make a dictionary of the expected samples
    found = dict()
    for sample in js["Samples"]:
        found[sample["Sample_Name"]] = 0

    # Get file list on the mounted webdav drive
    files = os.listdir(webdav)

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
            # print(sample, "OK")
            found[sample_name] = 1

    for sample, was_found in found.items():
        if was_found != 1:
            print("grep", sample, "basespace-copy-data.sh")

    # Now check if all uploaded files have a good file size on the webdav server
    print("============ CHECK FILE SIZE ===========")

    # Connect to all files and check the file size
    for myfile in files:
        sample, rest = myfile.split("_L001")
        link = args.webdav_url + args.run + "/data/" + myfile
        cmd = "curl -sI --netrc " + link 
        p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
        info = p.communicate()[0]
        info = info.split("\r\n")
	for entry in info:
            if entry.startswith("Content-Length"):
                fs = int(entry.replace("Content-Length: ", ""))
                break

        if fs < 1000:
            print("grep", myfile, "basespace-copy-data.sh")

