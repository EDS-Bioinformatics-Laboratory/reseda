import sys
import regex
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import IUPACData
import editdistance
from sequences import *

# Usage: python motif-search.py sequences.fasta
# Motifs are defined in this script
# Script prints BED entries to stdout

# Reads a fasta file and searches for (a set of) motifs

########## Config ##########

# BED format is 0-based
# GATK variant list might be 1-based
coordSystem = 0


######## Functions #########

########### Main ###########

motifs=["WRCY"]  # RGYW is reverse complement of WRCY
#motifs=["CGCCGCCAGCTCACC", "TTGCCCTCAACGACCACTTT", "CCTGCTTCTCCTCAGCTTCAG"]  # ACTB, GAPDH, HPRT1
# motifs=readMotifsFromFile("immuno-primers.csv")

# Check if an argument was given to this script
if len(sys.argv) < 2:
    sys.exit('Usage: %s sequences.fasta' % sys.argv[0])

# Add reverse complement of motifs to the list
motifs_comrev=list()
for motif in motifs:
    motifs_comrev.append(comrev(motif))
motifs = motifs + motifs_comrev

# file with sequences (e.g. test.fasta)
fastaFile = sys.argv[1]

# Translate motifs to regular expressions
myregex = dict()
for motif in motifs:
    myregex[motif] = motifToRegex(motif, 2)  # also adds in-exact matching pattern

# Open fasta file and search for all motifs in each sequence
for record in SeqIO.parse(open(fastaFile, "rU"), "fasta") :
    foundMatch = 0
    sys.stderr.write('Reading sequence record... ')
    sequence = str(record.seq).upper()
    sys.stderr.write('done\n')
    for motif, pat in myregex.iteritems():
        motifLength = len(motif)
        if motif in motifs_comrev:
            strand = "-"
        else:
            strand = "+"

        # Following code is from the method nt_search in BioPython (want to print result immediately iso putting it into a variable)
        # Replaced the re package with the regex package
        pos = -1
        l = len(sequence)
        while True:
            pos += 1
            s = sequence[pos:]
            m = regex.search(pat, s, regex.BESTMATCH)
            if not m:
                break
            pos += int(m.start(0))

            # Output BED entries: chr start end name score strand
            distance = editdistance.eval(motif, sequence[pos:pos+motifLength])
            # print record.id, pos+coordSystem, pos+motifLength+coordSystem, motif, sequence[pos:pos+motifLength], distance, strand, sequence[pos:]  # just space delimited
            print record.id, pos+coordSystem, pos+motifLength+coordSystem, sequence[pos:pos+motifLength], 1000, strand   # BED entry
            foundMatch = 1
            break                       # only record the best hit with this motif
    if foundMatch == 0:
        print record.id, "no-match"
