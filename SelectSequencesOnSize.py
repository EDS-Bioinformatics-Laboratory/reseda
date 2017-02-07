from __future__ import print_function
import sys

if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("Usage: " + sys.argv[0] + "min_length in.sam out.sam")
    (min_length, samIn, samOut) = sys.argv[1:]
    min_length = int(min_length)

    try:
        fhIn = open(samIn)
        fhOut = open(samOut, "w")
        fhRep = open(samIn + ".seqlength.report", "w")
    except:
        sys.exit("cannot read or write file: " + samIn + " " + samOut)

    total = 0
    kept = 0

    for line in fhIn:
        # if (m/^\@/) {print;} else { @c=split(/\s+/); print if length($c[9]) > 299; }
        line = line.rstrip()
        if line.startswith("@"):
            print(line, file=fhOut)
        else:
            total += 1
            c = line.split()
            if len(c[9]) >= min_length:
                kept += 1
                print(line, file=fhOut)

    print(kept, total, 100.0 * kept / total, file=fhRep)

    fhIn.close()
    fhOut.close()
    fhRep.close()
