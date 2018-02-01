'''
RESEDA - REPertoire SEquencing Data Analysis
Copyright (C) 2016 Barbera DC van Schaik

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from __future__ import print_function
import sys
import argparse
import gzip
import regex
from Bio import SeqIO
from sequences import *


def getVmotifs(cellType):
    '''
    Description: Retrieve motifs for V section
    In: cell type
    Out: list with V motifs
    '''
    motifs = list()

    if cellType == "IGH_HUMAN":
        refFile = "ref.table.heavy.csv"
    elif cellType == "TRB_HUMAN":
        refFile = "ref.table.TCRb.csv"
    elif cellType == "IGK_HUMAN":
        refFile = "ref.table.BCRk.csv"
    elif cellType == "IGL_HUMAN":
        refFile = "ref.table.BCRl.csv"
    elif cellType == "TRA_HUMAN":
        refFile = "ref.table.TCRa.csv"
    elif cellType == "IGH_MOUSE":
        refFile = "ref.table.mouse.heavy.csv"
    elif cellType == "IGK_MOUSE":
        refFile = "ref.table.mouse.BCRk.csv"
    elif cellType == "IGL_MOUSE":
        refFile = "ref.table.mouse.BCRl.csv"
    elif cellType == "TRB_MOUSE":    # Not tested on Miseq yet
        refFile = "ref.table.mouse.TCRb.csv"
    elif cellType == "TRA_MOUSE":    # Not tested yet and probably needs to be changed
        refFile = "ref.table.mouse.TCRa.csv"
    else:
        sys.exit("Cell type " + cellType + " is not implemented yet.\n" + usage)

    try:
        fhRef = open(refFile, "r")
    except:
        sys.exit("cannot open file:" + refFile)

    # Determine which column contains the sequence
    header = fhRef.readline()
    header = header.rstrip()
    header = header.replace("\"", "")
    cHeader = header.split(",")
    idx = -1
    for i in range(len(cHeader)):
        if cHeader[i] == "seq":
            idx = i

    if idx == -1:
        sys.exit("Couldn't find the 'seq' column in the reference file: " + refFile)

    # Read rest of the file
    for line in fhRef:
        line = line.rstrip()
        c = line.split(",")
        seq = c[idx]
        if cellType.startswith("TR"):   # TRB, TRA human
            motifs.append(seq[99:104])
        elif cellType.startswith("IGH"):  # IGH human
            motifs.append(seq[97:102])
        elif cellType == "IGL_HUMAN" or cellType == "IGK_HUMAN":
            motifs.append(seq[98:103])
        elif cellType == "IGL_MOUSE" or cellType == "IGK_MOUSE":
            motifs.append(seq[97:102])
        else:                           # Guess for new cell type
            motifs.append(seq[98:103])

    motifs = list(set(motifs))      # make list with motifs unique

    return(motifs)


def getJmotifs(cellType):
    '''
    Description: Return string with J motif
    In: cell type
    Out: list with J motifs
    '''
    if cellType == "IGH_HUMAN":
        return(["VTVS"])
    elif cellType == "IGH_MOUSE":
        return(["VTVS", "LTVS"])   # This needs to be verified by Sabrina and/or Giulia
    elif cellType == "TRB_HUMAN":
        return(["FG.G"])
    elif cellType == "TRB_MOUSE":
        return(["FG.G"])
    elif cellType == "IGK_HUMAN":
        return(["FG.G"])
    elif cellType == "IGK_MOUSE":
        return(["FG.G"])
    elif cellType == "IGL_HUMAN":
        return(["FG.G"])
    elif cellType == "IGL_MOUSE":
        return(["FG.G"])
    elif cellType == "TRA_HUMAN":
        return(["FG.G", "FARG", "WGAG", "WGLG"])
    elif cellType == "TRA_MOUSE":
        return(["FG.G", "FARG", "WGAG", "WGLG"])
    else:
        sys.exit("Cell type " + cellType + " is not implemented yet.\n" + usage)


def getMotifs(cellType, mismatches):
    '''
    Description: retrieve V and J motifs and concatenate them
    In: string cellType, int mismatches (usually 0 or 1)
    Out: string combinedMotifs (a regular expression)
    '''

    if mismatches == 0 and type(mismatches) == type(10):
        mismatches = ""
    elif mismatches > 0 and type(mismatches) == type(10):
        mismatches = "{e<=" + str(mismatches) + "}"
    else:
        raise TypeError('wrong input for mismatches')

    motifsV = getVmotifs(cellType)
    motifsJ = getJmotifs(cellType)

    # print("V:", ",".join(motifsV))
    # print("J:", ",".join(motifsJ))

    motifs = list()
    for v in motifsV:

        for j in motifsJ:
            if v == "" or j == "":  # skip empty values
                continue

            # v = v[:-1] + "(" + v[-1]    # put a "(" between second last and last character
            motifs.append(v + ".+?" + j)

    combinedMotifs = ".+(" + "|".join(motifs) + ")" + mismatches
    # print(combinedMotifs)

    return(combinedMotifs)


def getAlternativeVmotifs(cellType, mismatches):
    '''
    Description: retrieve V motifs and add a wildcard to it (10 times .)
    In: string cellType, int mismatches (usually 0 or 1)
    Out: string combinedMotifs (a regular expression)
    '''

    if mismatches == 0 and type(mismatches) == type(10):
        mismatches = ""
    elif mismatches > 0 and type(mismatches) == type(10):
        mismatches = "{e<=" + str(mismatches) + "}"
    else:
        raise TypeError('wrong input for mismatches')

    motifsV = getVmotifs(cellType)

    motifs = list()
    for v in motifsV:
        if v == "":  # skip empty values
            continue

        motifs.append(v + ".{13}")

    combinedMotifs = ".+(" + "|".join(motifs) + ")" + mismatches

    return(combinedMotifs)


def getAlternativeJmotifs(cellType, mismatches):
    '''
    Description: retrieve J motifs and add a wildcard in front of it (10 times .)
    In: string cellType, int mismatches (usually 0 or 1)
    Out: string combinedMotifs (a regular expression)
    '''

    if mismatches == 0 and type(mismatches) == type(10):
        mismatches = ""
    elif mismatches > 0 and type(mismatches) == type(10):
        mismatches = "{e<=" + str(mismatches) + "}"
    else:
        raise TypeError('wrong input for mismatches')

    motifsJ = getJmotifs(cellType)

    motifs = list()
    for j in motifsJ:
        if j == "":  # skip empty values
            continue

        motifs.append(".{14}" + j)

    combinedMotifs = ".+(" + "|".join(motifs) + ")" + mismatches

    return(combinedMotifs)


def getImgtMotifs():
    '''
    Description: Return regular expression for Cys-Phe/Trp motif or Cys-Val
    In: string cellType
    Out: motif
    '''

    # Depending on cellType CDR3 ends with F, W or V
    combinedMotifs = ".+(C.+?[FWV]..)"
    # combinedMotifs = ".+(C.+[FWV]..)"

    return(combinedMotifs)


def extractCDR3(cellType, peptide, p):
    '''
    Description: extract the CDR3 from a peptide sequence
    In: peptide sequence (string), regular expression (regex.compile object) V+J, regular expression for just the V (regex.compile object)
    Out: CDR3 sequence (string), aa_pos: list with start (int) and end (int) positions in the peptide sequence
    '''

    cdr3pep = None
    aa_pos = list()

    m = p.search(str(peptide))

    if m is not None:   # a match is found
        # Extract CDR3 peptide sequence
        cdr3pep = m.group(1)
        aa_pos = list(m.span(1))

        # First check if there is a Cys in the peptide, in that case report CDR3 from there
        if "C" in cdr3pep:
            c_pos = cdr3pep.find("C")
            cdr3pep = cdr3pep[c_pos:]
            aa_pos[0] = aa_pos[0] + c_pos
        else:
            # Remove the first 6, 5 or 4 aminoacids of the CDR3
            if cellType.startswith("IGH"):  # IGH
                cdr3pep = cdr3pep[6:]
                aa_pos[0] = aa_pos[0] + 6
            elif cellType.startswith("IG"):  # IGL, IGK
                cdr3pep = cdr3pep[6:]
                aa_pos[0] = aa_pos[0] + 6
            else:                           # TRB, TRA
                cdr3pep = cdr3pep[4:]
                aa_pos[0] = aa_pos[0] + 4

        # Remove the last 2 aminoacids of the CDR3
        cdr3pep = cdr3pep[:-2]
        aa_pos[1] = aa_pos[1] - 2

    return(cdr3pep, aa_pos)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Translates nucleotide to protein and extracts CDR3')
    parser.add_argument('-c', '--celltype', default='IGH_HUMAN', type=str, help='Cell type: (IGH|TRB|IGK|IGL|TRA)_HUMAN (default: %(default)s)')
    parser.add_argument('-m', '--mismatches', default=0, type=int, help='Allowed nr of mismatches in motif search (default: %(default)s)')
    parser.add_argument("fastq_files", type=str, nargs='+', help='Path(s) to fastq file(s)')
    args = parser.parse_args()

    cellType = args.celltype

    # Get all the motifs to search for V .* J, define mismatches (usually 0 or 1)
    motif = getMotifs(cellType, args.mismatches)
    p = regex.compile(motif, regex.BESTMATCH)

    # Alternative motif: search for V motif plus 10 amino acids and 10aa + J motif
    motifAltV = getAlternativeVmotifs(cellType, args.mismatches)
    p_alt_v = regex.compile(motifAltV, regex.BESTMATCH)
    motifAltJ = getAlternativeJmotifs(cellType, args.mismatches)
    p_alt_j = regex.compile(motifAltJ, regex.BESTMATCH)

    # Check for an extra motif Asn-X-Ser/Thr (X is not Proline)
    p_extra = regex.compile("N[^P][ST]")

    # Pattern for stop codon and untranslated codons
    p_stop = regex.compile("\*")
    p_x = regex.compile("X")

    # Open fastq file(s) and search for patterns
    for inFile in args.fastq_files:
        outFile = inFile + "-" + cellType + "-CDR3.csv"
        altVFile = inFile + "-" + cellType + "-alt-V-CDR3.csv"
        altJFile = inFile + "-" + cellType + "-alt-J-CDR3.csv"
        extraFile = inFile + "-" + cellType + "-extra.txt"
        rawFile = inFile + "-" + cellType + ".csv"
        repFile = inFile + "-" + cellType + "-report.txt"
        stopFile = inFile + "-" + cellType + "-discarded-stop-codon.txt"
        uncalledFile = inFile + "-" + cellType + "-discarded-uncalled-bases.txt"
        nocdr3File = inFile + "-" + cellType + "-discarded-no-cdr3.txt"
        try:
            fhIn = gzip.open(inFile, "rb")
        except:
            sys.exit("cannot open file: " + inFile)
        try:
            fhOut = open(outFile, "w")
        except:
            sys.exit("cannot write to file: " + outFile)
        try:
            fhRaw = open(rawFile, "w")
        except:
            sys.exit("cannot write to file: " + rawFile)
        try:
            fhRep = open(repFile, "w")
        except:
            sys.exit("cannot write to file: " + repFile)
        try:
            fhStop = open(stopFile, "w")
        except:
            sys.exit("cannot write to file: " + stopFile)
        try:
            fhUncalled = open(uncalledFile, "w")
        except:
            sys.exit("cannot write to file: " + uncalledFile)
        try:
            fhNoCdr3 = open(nocdr3File, "w")
        except:
            sys.exit("cannot write to file: " + nocdr3File)
        try:
            fhExtra = open(extraFile, "w")
        except:
            sys.exit("cannot write to file:" + extraFile)
        try:
            fhAltV = open(altVFile, "w")
        except:
            sys.exit("cannot write to file:" + altVFile)
        try:
            fhAltJ = open(altJFile, "w")
        except:
            sys.exit("cannot write to file:" + altJFile)

        # Container to count stuff
        count_stuff = dict()
        count_accs_alt_v = list()
        count_accs_alt_j = list()

        for record in SeqIO.parse(fhIn, "fastq"):

            count_stuff["1. Total reads"] = count_stuff.get("1. Total reads", 0) + 1

            # Translate sequence to protein
            translations = nucToPeptide(str(record.seq))

            # Was a CDR3 found?
            cdr3_found = False

            # Find motif (V, J and everything in between)
            # TO DO: sometimes multiple patterns can match which results in a much longer CDR3. Needs a fix
            for i in range(len(translations)):
                (cdr3pep, aa_pos) = extractCDR3(cellType, str(translations[i]), p)

                # Search for motif in the protein translation
                for m_extra in p_extra.finditer(str(translations[i])):
                    print(record.id, str(i), m_extra.group(0), m_extra.span(), file=fhExtra)

                if cdr3pep is not None:
                    cdr3_found = True

                    # Extract CDR3 nucleotide sequence
                    nt_start = aa_pos[0] * 3
                    nt_end = aa_pos[1] * 3
                    if i < 3:                   # retrieve nucleotide reading frame +
                        tmp_seq = record.seq[i:]
                    else:                       # retrieve nucleotide reading frame of reverse complement
                        tmp_seq = comrev(record.seq)[i - 3:]
                    cdr3nuc = tmp_seq[nt_start:nt_end]

                    # Convert fastq quality scores to phred scores
                    quality_scores = record.letter_annotations["phred_quality"]
                    if i < 3:
                        tmp_qual = quality_scores[i:]
                    else:
                        tmp_qual = list(reversed(quality_scores))[i - 3:]
                    cdr3_quality_scores = tmp_qual[nt_start:nt_end]

                    # Basic stats CDR3 quality
                    if len(cdr3_quality_scores) > 0:
                        cdr3_qual_min = min(cdr3_quality_scores)
                        cdr3_qual_max = max(cdr3_quality_scores)
                        cdr3_qual_avg = round(sum(cdr3_quality_scores) / float(len(cdr3_quality_scores)), 1)
                    else:
                        cdr3_qual_min = 0
                        cdr3_qual_max = 0
                        cdr3_qual_avg = 0
                        cdr3_quality_scores = "-"

                    # Convert the quality scores to a string
                    quality_scores = str(quality_scores)
                    quality_scores = quality_scores.replace(",", "").replace("[", "").replace("]", "")
                    cdr3_quality_scores = str(cdr3_quality_scores)
                    cdr3_quality_scores = cdr3_quality_scores.replace(",", "").replace("[", "").replace("]", "")

                    # Check for stop codons and untranslated codons
                    if p_stop.search(cdr3pep) is not None:    # stop codon in cdr3peptide
                        count_stuff["2. Discarded reads stop codon in CDR3"] = count_stuff.get("2. Discarded reads stop codon in CDR3", 0) + 1
                        print("\t".join([record.id, str(i), str(record.seq), str(translations[i]), quality_scores]), file=fhStop)
                    elif p_x.search(cdr3pep) is not None:     # check for uncalled bases
                        count_stuff["3. Discarded reads uncalled bases in CDR3"] = count_stuff.get("3. Discarded reads uncalled bases in CDR3", 0) + 1
                        print("\t".join([record.id, str(i), str(record.seq), str(translations[i]), quality_scores]), file=fhUncalled)
                    else:                                 # Correct CDR3 peptide found without uncalled bases or stop codons
                        print("\t".join([record.id, str(i), str(cdr3pep), str(cdr3nuc), str(cdr3_qual_min), str(cdr3_qual_max), str(cdr3_qual_avg), cdr3_quality_scores, str(nt_start), str(nt_end), str(len(record.seq))]), file=fhOut)
                        print("\t".join([record.id, str(i), str(record.seq), str(translations[i]), quality_scores]), file=fhRaw)

                        count_stuff["4. Reads with CDR3"] = count_stuff.get("4. Reads with CDR3", 0) + 1
                        break

            # Nothing found: try searching for V + 10 amino acids
            if cdr3_found is False:
                for i in range(len(translations)):
                    (cdr3pep, aa_pos) = extractCDR3(cellType, str(translations[i]), p_alt_v)
                    if cdr3pep is not None:
                        cdr3_found = True
                        count_accs_alt_v.append(record.id)
                        print("\t".join([record.id, str(i), str(record.seq), str(cdr3pep)]), file=fhAltV)

            # Still nothing found: try searching for 10 amino acids + J
            if cdr3_found is False:
                # Make an alternative clones file
                for i in range(len(translations)):
                    (cdr3pep, aa_pos) = extractCDR3(cellType, str(translations[i]), p_alt_j)
                    if cdr3pep is not None:
                        cdr3_found = True
                        count_accs_alt_j.append(record.id)
                        print("\t".join([record.id, str(i), str(record.seq), str(cdr3pep)]), file=fhAltJ)

            # Could not find anything
            if cdr3_found is False:
                print("\t".join([record.id, str(i), str(record.seq), str(translations)]), file=fhNoCdr3)

        # Make report
        print("Motifs:", motif, file=fhRep)

        count_stuff["5. Reads with only V motif"] = len(set(count_accs_alt_v))
        count_stuff["6. Reads with only J motif"] = len(set(count_accs_alt_j))

        total = count_stuff["1. Total reads"]
        for key, value in sorted(count_stuff.iteritems()):
            perc = 100.00 * value / total
            print("\t".join([key, str(value), str(perc) + "%"]), file=fhRep)

        fhIn.close()
        fhOut.close()
        fhRaw.close()
        fhRep.close()
        fhStop.close()
        fhUncalled.close()
        fhNoCdr3.close()
        fhExtra.close()
        fhAltV.close()
        fhAltJ.close()
