from __future__ import print_function
import os.path

if __name__ == '__main__':
    report_file = "run17-report-all.csv"
    file_list = "ls-assembled.txt"

    # Which MID is the correct one for a sample?
    fh = open(report_file)
    header = fh.readline()
    correct_mid = dict()
    for line in fh:
        line = line.rstrip()
        line = line.split(" ")
        correct_mid[line[0]] = line[4]
    fh.close()

    # Select file paths with requested sample and the correct MID
    fh = open(file_list)
    for file_path in fh:
        file_path = file_path.rstrip()
        basename = os.path.basename(file_path)
        dirname = os.path.dirname(file_path)
        sample, mid = basename.split("_L001.assembled-")
        mid = mid.replace(".fastq.gz", "")
        mid = mid.replace("_fastqc.zip", "")
        try:  # sample should be in the dictionary, otherwise skip this line
            get_mid = correct_mid[sample]
        except:
            continue
        if mid == get_mid:    # if mid is the correct mid, then print the path
            print("mv -v", file_path, dirname + "/fastqc-correct-mid/")
    fh.close()
