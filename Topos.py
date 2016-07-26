from __future__ import print_function
import sys
import os
import subprocess

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

def uploadToken (poolname, f):
    '''
    Description: upload file as token to topos
    In: filepath (str)
    Out: token number (str)
    '''
    cmd = "topos uploadFileAsToken " + poolname + " " + f
    syscall = os.popen(cmd)
    token = syscall.read()
    syscall.close()
    return(token)

if __name__ == '__main__':
    if len(sys.argv)<3:
        sys.exit("Usage: " + sys.argv[0] + " poolname files-to-upload")

    poolname = sys.argv[1]
    for myfile in sys.argv[2:]:
        token = uploadToken(poolname, myfile)
        print(token)