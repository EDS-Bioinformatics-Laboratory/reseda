from __future__ import print_function
import os
import argparse

'''
Description: get file names of all fastq files in the basespace directory
In: fill in the sampledir and beehub url below
Out: curl commands on stdout
'''

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Print curl commands to transfer Basespace data to Webdav server')
    parser.add_argument('-r', '--run', default='runNN-yyyymmdd-miseq', type=str, help='Run directory name (default: %(default)s)')
    parser.add_argument('-b', '--basespace_dir', default='/data/home/bioinfo/basespace/Projects/', type=str, help='Basespace directory  (default: %(default)s)')
    parser.add_argument("sub_dirs", type=str, nargs='+', help='Basespace sub-directorie(s)')
    args = parser.parse_args()

    if args.run == 'runNN-yyyymmdd-miseq':
        parser.print_help()
        exit()

    myurl = "https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/" + args.run + "/data/"

    fhOut = open("basespace-copy-data.sh", "w")
    for mydir in args.sub_dirs:
        # syscall = os.popen("ls " + args.basespace_dir + mydir + "/*/*.fastq.gz")
        # syscall = os.popen("find " + args.basespace_dir + mydir + " -maxdepth 5 -mtime 5 -regex '.*.fastq.gz'")
        syscall = os.popen("ls " + args.basespace_dir + mydir + "/Samples/*/Files/*.fastq.gz")

        for line in syscall:
            line = line.rstrip()
            print('curl -T "' + line + '" --netrc', myurl, "; wait", file=fhOut)
    fhOut.close()
    print("Wrote basespace-copy-data.sh to disk")
