from __future__ import print_function
import sys


def fixCoordinates(coord):
    '''
    Description: merge coordinates that are next to each other
    In: list of tuples
    Out: list of tuples
    '''
    newcoord = list()
    (prev_start, prev_end) = coord[0]
    for (start, end) in coord[1:]:
        if start == prev_end:
            # merge with previous
            prev_end = end
        else:
            # store coordinates in new list
            newcoord.append((prev_start, prev_end))
            (prev_start, prev_end) = (start, end)
    newcoord.append((prev_start, prev_end))

    return(newcoord)


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
            if character == "H" or character == "D":   # do not use hard clipped numbers or deletions
                number = ""
                continue

            end = start + number
            if character == "M" or character == "I":
                coord.append((start, end))
            start = end

            # reset number
            number = ""

    # Merge coordinates if they are next to each other
    coord = fixCoordinates(coord)

    return(coord)


def parseSam(f):
    '''
    Description: gets number of mutations in the aligned part of the query sequence
    In: name of sam file
    Out: name of file with mutations
    '''
    try:
        fh = open(f)
        outfile = f + '.mut.txt'
        fhOut = open(outfile, "w")
    except:
        sys.exit("cannot open or write file: " + f)

    print("acc cigar start.pos end.pos mut.count mut.perc align.length align.seq", file=fhOut)

    for line in fh:
        # Skip header
        if line.startswith("@"):
            continue

        # Read alignment
        line = line.rstrip()
        line = line.split()
        try:
            cigar = parseCigar(line[5])
        except:
            continue
        acc = line[0]
        seq = line[9]
        if len(cigar) > 1:
            print("WARNING: multiple parts aligned:", acc, line[5], cigar, "\n", seq, "\n-----------------------------------------------------------------------")
        for (start, end) in cigar:
            subseq = seq[start:end]
            countEqual = subseq.count("=")
            countMut = len(subseq) - countEqual
            try:
                percMut = 100 * countMut / len(subseq)
            except:
                print("WARNING: aligned part has zero length:", acc, line[5], seq)
                exit()
            print(acc, line[5], start, end, countMut, percMut, len(subseq), subseq, file=fhOut)

    print("Wrote", outfile, "to disk")
    return(outfile)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage:", sys.argv[0], "*.sam")
        exit()

    for samfile in sys.argv[1:]:
        parseSam(samfile)
