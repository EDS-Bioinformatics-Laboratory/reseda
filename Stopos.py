from __future__ import print_function
import os
import subprocess


def deleteToken(poolname, token):
    '''
    Description: delete token from pool
    In: poolname (str), token (str)
    Out: -
    '''
    cmd = ["stoposclient", "remove", "-p", poolname, token]
    rc = subprocess.call(cmd)
    if os.environ.get('STOPOS_RC') != "OK":
        print("WARNING: couldn't delete token from pool:", poolname, token)


def uploadToken(poolname, f):
    '''
    Description: upload file as token to topos
    In: filepath (str)
    Out: token number (str)
    '''
    cmd = ["stoposclient", "add", "-p", poolname, f]
    rc = subprocess.call(cmd)
    if os.environ.get('STOPOS_RC') != "OK":
        print("WARNING: couldn't add token to pool:", poolname, f)
    token = os.environ.get("STOPOS_KEY")
    return(token)


if __name__ == '__main__':
    print("Nothing in main")
