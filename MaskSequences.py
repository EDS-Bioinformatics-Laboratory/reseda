import sys
import gzip
from MutationsFromSam import parseCigar


def maskSequence(seq, coords, letter="n"):
    '''
    Description: mask a sequence
    In: seq (str), coords (list of tuples with start and end), letter (e.g. "N")
    Out: maskedseq (str)
    '''
    for start, end in coords:
        upstream = seq[:start]
        subseq = seq[start:end]
        downstream = seq[end:]
        mask = letter * len(subseq)
        seq = upstream + mask + downstream

    return(seq)


def maskSam(samfile):
    '''
    Description: masks sequences using the cigar string
    In: sam file (filename, str)
    Out: fastq file (filename, str)
    '''
    maskedsam = samfile.replace(".sam", ".masked.fastq.gz")
    try:
        fhIn = open(samfile)
        # fhOut = gzip.open(maskedsam, "w")
    except:
        # sys.exit("Cannot open or write file " + samfile + " " + maskedsam)
        sys.exit("Cannot open file " + samfile)

    with gzip.open(maskedsam, "wt") as fhOut:
        for line in fhIn:
            line = line.rstrip()
            if line.startswith("@"):
                # print(line, file=fhOut)
                pass
            else:
                line = line.split()
                acc = line[0]
                cigar = line[5]
                seq = line[9]
                qual = line[10]
                try:
                    seq = maskSequence(seq, parseCigar(cigar))
                except:
                    pass
                print("@" + acc, file=fhOut)
                print(seq, file=fhOut)
                print("+", file=fhOut)
                print(qual, file=fhOut)
                if(len(seq) != len(qual)):
                    print(acc, "length not equal")

    fhIn.close()
    fhOut.close()
    return(maskedsam)


def test_mask():
    seq = "abcdefg"
    coords = [(2,5)]  # "cde"
    print("Test one okay?", maskSequence(seq, coords, "N") == 'abNNNfg')
    coords = [(1,3),(4,6)]  # "bc" and "ef"
    print("Test two okay?", maskSequence(seq, coords, "N") == 'aNNdNNg')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage:", sys.argv[0], "*.sam")
        print("Example: python MaskSequences.py *.sam")
        print("Will create new files with the extension .masked.sam")
        exit()

    for samfile in sys.argv[1:]:
        maskedsam = maskSam(samfile)
        print("Masked sam file:", maskedsam)
