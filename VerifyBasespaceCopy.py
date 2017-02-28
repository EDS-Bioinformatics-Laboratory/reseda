from __future__ import print_function
import os

if __name__ == '__main__':
    mnt_dir = "/mnt/immunogenomics/RUNS//run13-20170224-miseq/data"
    files = os.listdir(mnt_dir)

    check = dict()
    for myfile in files:
        sample, rest = myfile.split("_L001")
        check[sample] = check.get(sample, list())
        check[sample].append(myfile)

    for sample, pairs in check.items():
        if len(pairs) < 2:
            if "R1" in pairs[0]:
                print("grep", sample, "basespace-copy-data.sh | grep R2")
            elif "R2" in pairs[0]:
                print("grep", sample, "basespace-copy-data.sh | grep R1")
