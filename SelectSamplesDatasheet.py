from __future__ import print_function
import sys
import re

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " sample-list.csv datasheet.csv")

    sample_list = sys.argv[1]
    data_sheet = sys.argv[2]

    # Open file with requested samples and put them in a list

    try:
        fh = open(sample_list)
    except:
        sys.exit("cannot open file: " + sample_list)

    samples = list()
    header = fh.readline()
    header = header.rstrip()
    header = header.split(",")
    for line in fh:
        line = line.rstrip()
        line = line.split(",")
        sample = line[header.index("SampleID")]
        sample = re.sub("_S\d+$", "", sample)
        samples.append(sample)

    fh.close()
    # print(samples)

    # Open datasheet and keep samples that are in the list from above
    try:
        fh = open(data_sheet)
    except:
        sys.exit("cannot open file: " + data_sheet)

    data = 0
    for line in fh:
        line = line.rstrip()
        if data == 0:                   # Meta data before the datasets
            print(line)
            if line.startswith("[Data]"):
                data = 1
                header = fh.readline()
                header = header.rstrip()
                print(header)
                header = header.split(",")
        else:                           # Sample information
            c = line.split(",")
            sample = c[header.index("Sample_Name")]
            if sample in samples:       # Only print samples that are in the samples list
                print(line)
    fh.close()
