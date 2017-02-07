from __future__ import print_function
import sys
from Topos import *

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit("Usage: " + sys.argv[0] + " poolname files-to-upload")

    poolname = sys.argv[1]
    for myfile in sys.argv[2:]:
        token = uploadToken(poolname, myfile)
        print(token)
