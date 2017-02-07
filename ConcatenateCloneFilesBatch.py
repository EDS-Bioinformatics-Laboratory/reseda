from __future__ import print_function
# import subprocess


def executeCmd(cmd):
    '''
    Description: execute a command and check if everything went well
    In: command (str)
    Out: - an exception will be raised when command was not successful
    '''
    print(cmd)
    # cmd = cmd.split()
    # rc = subprocess.call(cmd)
    # if rc != 0:
    #     print("something went wrong with: " + " ".join(cmd))


if __name__ == '__main__':
    mydir = "/mnt/immunogenomics/RUNS/run12-20170127-miseq/results-tbcell/final/correct-mid/"
    celltypes = ["IGH_HUMAN", "IGH_MOUSE", "IGK_HUMAN", "IGK_MOUSE", "IGL_HUMAN", "IGL_MOUSE", "TRA_HUMAN", "TRA_MOUSE", "TRB_HUMAN", "TRB_MOUSE"]

    for celltype in celltypes:
        cmd = "python ConcatenateCloneFiles.py " + mydir + "*" + celltype + "-clones-subs.csv"
        executeCmd(cmd)
        cmd = "wait"
        executeCmd(cmd)
        cmd = "mv run-clones_subs.csv run-clones_subs-" + celltype + ".csv"
        executeCmd(cmd)
        cmd = "python ConcatenateCloneFiles.py " + mydir + "*" + celltype + "*.rr.clones_subs.csv"
        executeCmd(cmd)
        cmd = "wait"
        executeCmd(cmd)
        cmd = "mv run-clones_subs.csv run-clones_subs-" + celltype + "-after-reassignment.csv"
        executeCmd(cmd)
