from __future__ import print_function
import sys


def parseCigar(cigarstring):
    '''
    Description: split cigarstring and store in list of tuples [(number, letter)] and convert it to coordinates in sequence
    '''
    c = list()
    coord = list()
    number = ""
    start = 0
    end = 0
    for character in cigarstring:
        try:
            x = int(character)
            number += character
        except:
            number = int(number)
            c.append((number, character))

            # convert to coordinates if part of the alignment
            end = start + number
            if character == "M":
                coord.append((start, end))
            start = end

            # reset number
            number = ""

    return(coord)


def parseSam(f):
    '''
    Description: gets number of mutations in the aligned part of the query sequence
    In: name of sam file
    Out: -
    '''
    try:
        fh = open(f)
    except:
        sys.exit("cannot open file: " + f)

    for line in fh:
        # Skip header
        if line.startswith("@"):
            continue

        # Read alignment
        line = line.rstrip()
        line = line.split()
        cigar = parseCigar(line[5])
        seq = line[9]
        print(line[5], cigar, seq)
        for (start, end) in cigar:
            print("  " + seq[start:end])


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage:", sys.argv[0], "*.sam")
        exit()

    for samfile in sys.argv[1:]:
        parseSam(samfile)
