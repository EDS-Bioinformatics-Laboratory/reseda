testmode = 0    # is the script in testmode? 0=no 1=yes

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import IUPACData
from Bio.Alphabet import generic_dna

# Required: directory "hg19" with fasta files per chromosome (from UCSC)

def readGenome ():
    """
    Read all chromosomes and put sequence into dictionary "genome"
    """
    genome = dict()
    chromosomes = range(1,23) + ["X", "Y", "M"]
    #chromosomes = [1, 3]
    for chrom in chromosomes:
        fastaFile = "hg19/chr" + str(chrom) + ".fa"
        for record in SeqIO.parse(open(fastaFile, "rU"), "fasta") :
            genome[record.id] = str(record.seq).upper()

    return genome

def motifToRegex (motif, mismatches):
    """
    Converts IUPAC motif to a regular expression (copied from nt_search method in BioPython)
    """
    pattern = ''
    for nt in motif:
        value = IUPACData.ambiguous_dna_values[nt]
        if len(value) == 1:
            pattern += value
        else:
            pattern += '[%s]' % value
    pattern = "(" + pattern + "){e<=" + str(mismatches) + "}"  # allow for mismatches
    return pattern

def complement(s): 
    """
    Return the complement nucleotides
    """
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y':'R', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)

def comrev(s):
    """
    Return reverse complement
    """
    return complement(s[::-1])


def readMotifsFromFile (motifFile):
    """
    Read a tab-delimited files with motifs to search for. Return a list with motifs
    """
    try:
        fh = open(motifFile)
    except:
        sys.exit("cannot open file", fileList)

    motifs = list()
    line = fh.readline() # skip header
    for line in fh:
        line = line.rstrip()
        c = line.split("\t")
        motifs.append(c[1].upper())
    return motifs


def nucToPeptide (seq):
    '''
    Description: translates a nucleotide sequence in all 6 reading frames
    In: string sequence
    Out: list with translations
    '''
    translations = list()
    frame = list()

    # Make sequence complement reverse
    seq = seq.upper()
    comrevSeq = comrev(seq)

    # Get all reading frames
    frame.append(Seq(seq, generic_dna))
    frame.append(Seq(seq[1:], generic_dna))
    frame.append(Seq(seq[2:], generic_dna))
    frame.append(Seq(comrevSeq, generic_dna))
    frame.append(Seq(comrevSeq[1:], generic_dna))
    frame.append(Seq(comrevSeq[2:], generic_dna))

    # Translate sequence to protein and find motif
    for i in range(len(frame)):
        translations.append(frame[i].translate())

    return(translations)


############### Tests ##############
if testmode==1:
    genome = readGenome()
    print genome["chr1"][0:50]
    