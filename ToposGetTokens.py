from __future__ import print_function
import sys
import os
import subprocess
import json

def getJob (poolname):
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

def deleteToken (poolname, token):
    '''
    Description: delete token from pool
    In: poolname (str), token (str)
    Out: -
    '''
    cmd = ["topos", "deleteToken", poolname, token]
    rc = subprocess.call(cmd)
    if rc != 0:
        print("WARNING: couldn't delete token from pool:", poolname, token)

def runJob (token, js):
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

def cleanUp ():
    '''
    Description: Clean up input and result files
    In: -
    Out: -
    '''
    cmd = ["rm", "-f", "*L001*"]
    print(cmd)
    rc = subprocess.call(cmd)

    cmd = ["rm", "-rf", "split/*L001*"]
    print(cmd)
    rc = subprocess.call(cmd)

    cmd = ["rm", "-f", "final/*L001*"]
    print(cmd)
    rc = subprocess.call(cmd)

    cmd = ["rm", "-f", "final/correct-mid/*L001*"]
    print(cmd)
    rc = subprocess.call(cmd)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("Usage: " + sys.argv[0] + " POOLNAME")

    poolname = sys.argv[1]

    # Prepare and run jobs. Script dies when there are no jobs left.
    while True:
        token,js = getJob(poolname)
        rc = runJob(token, js)

        # delete token from pool if run was successful
        if rc == 0:
            deleteToken(poolname, token)
        else:
            print("ERROR: job was not successful", token, js)

        # clean up result files before the next job
        cleanUp()
