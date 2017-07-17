from __future__ import print_function
import sys
import os
import shutil
import subprocess
import json
from Topos import *


def getJob(poolname):
    '''
    Description: retrieve next job from topos
    In: poolname
    Out: token (str), js (json)
    '''

    # Get token number
    cmd = "topos nextTokenWithLock " + poolname + " 3600"
    syscall = os.popen(cmd)
    token = syscall.readline()
    syscall.close()

    # Exit when there are no tokens left
    if token == "":
        print("FINISHED")
        exit()

    print("token", token)

    # Get token content
    cmd = "topos getToken " + poolname + " " + token
    syscall = os.popen(cmd)
    text = syscall.read()
    js = json.loads(text)
    syscall.close()

    return(token, js)


def runJob(token, js):
    '''
    Description: Prepare and run job
    In: token (str), js (json)
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
    cmd = ["./execute-all.sh", js["run"], js["mids"], js["organism"], js["cell"], js["celltype"]]
    print(cmd)
    rc = subprocess.call(cmd)
    return(rc)


def cleanUp():
    '''
    Description: Clean up input and result files
    In: -
    Out: -
    '''

    myfiles = os.listdir(".")
    for myfile in myfiles:
        if "L001" in myfile:
            if os.path.isfile(myfile):
                os.remove(myfile)
            else:
                shutil.rmtree(myfile, ignore_errors=True)

    myfiles = os.listdir("./split")
    for myfile in myfiles:
        if "L001" in myfile:
            shutil.rmtree("split/" + myfile, ignore_errors=True)

    myfiles = os.listdir("./split")
    for myfile in myfiles:
        if "L001" in myfile:
            os.remove("split/" + myfile)

    myfiles = os.listdir("./final")
    for myfile in myfiles:
        if "L001" in myfile:
            os.remove("final/" + myfile)

    myfiles = os.listdir("./final/correct-mid")
    for myfile in myfiles:
        if "L001" in myfile:
            os.remove("final/correct-mid/" + myfile)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " POOLNAME")

    poolname = sys.argv[1]

    # Prepare and run jobs. Script dies when there are no jobs left.
    while True:
        token, js = getJob(poolname)
        rc = runJob(token, js)

        # delete token from pool if run was successful
        if rc == 0:
            deleteToken(poolname, token)
        else:
            print("ERROR: job was not successful", token, js)

        # clean up result files before the next job
        cleanUp()

    cleanUp()
