from __future__ import print_function
import sys

if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("Usage: " + sys.argv[0] + " text-to-replace new-text token.json(s)")

    text_from = sys.argv[1]
    text_to = sys.argv[2]

    for myfile in sys.argv[3:]:
        try:
            fhIn = open(myfile)
            outfile = myfile + ".new.json"
            fhOut = open(outfile, "w")
        except:
            sys.exit("cannot open or write file: ", myfile)

        for line in fhIn:
            line = line.rstrip()
            line = line.replace(text_from, text_to)
            print(line, file=fhOut)

        fhIn.close()
        fhOut.close()
