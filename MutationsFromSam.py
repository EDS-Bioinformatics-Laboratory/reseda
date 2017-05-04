from __future__ import print_function
import sys


def parseSam(f):
    '''
    Description: gets number of mutations in the aligned part of the query sequence
    In: name of sam file
    Out: -
    '''
    print(f)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage:", sys.argv[0], "*.sam")
        exit()

    for samfile in sys.argv[1:]:
        parseSam(samfile)
