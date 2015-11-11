from __future__ import print_function
import sys
import gzip
import regex
from Bio import SeqIO
from sequences import *

usage = "Usage: " + sys.argv[0] + " (IGH|TRB|IGK|IGL|TRA)_HUMAN fastq-file(s)"

########### Functions ############

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
        if cellType.startswith("TR"):   # TRB, TRA
            motifs.append(seq[99:104])
        elif cellType.startswith("IGH"): # IGH
            motifs.append(seq[97:102])
        else:                           # IGK, IGL
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
    elif cellType == "TRB_HUMAN":
        return(["FG.G"])
    elif cellType == "IGK_HUMAN":
        return(["FG.G"])
    elif cellType == "IGL_HUMAN":
        return(["FG.G"])
    elif cellType == "TRA_HUMAN":
        return(["FG.G", "FARG", "WGAG", "WGLG"])
    else:
        sys.exit("Cell type " + cellType + " is not implemented yet.\n" + usage)

def getMotifs(cellType):
    '''
    Description: retrieve V and J motifs and concatenate them
    In: string cellType
    Out: list motifs
    '''
    motifsV = getVmotifs(cellType)
    motifsJ = getJmotifs(cellType)

    motifs = list()
    for v in motifsV:

        for j in motifsJ:
            if v == "" or j == "":  # skip empty values
                continue

            #v = v[:-1] + "(" + v[-1]    # put a "(" between second last and last character
            motifs.append(v + ".+?" + j)

    if cellType.startswith("TR"):                       # perform an exact match for TRA and TRB
        combinedMotifs = ".+(" + "|".join(motifs) + ")"
    else:                                               # allow max one mismatch for IGH, IGK, IGL
        combinedMotifs = ".+(" + "|".join(motifs) + "){e<=1}"
    print(combinedMotifs)

    return(combinedMotifs)

############# Main ###############

if len(sys.argv) < 3:
    sys.exit(usage)

cellType = sys.argv[1]

# Get all the motifs to search for V .* J
motif = getMotifs(cellType)

# Transform motif to regular expressions
p = regex.compile(motif, regex.BESTMATCH)

# Get V motifs (this is used to check if more than one motif is present)
p_v = regex.compile("(" + "|".join(getVmotifs(cellType)) + "){e<=0}", regex.BESTMATCH)

# Check for an extra motif
#p_extra = regex.compile("A[^P][ST]")

# Pattern for stop codon and untranslated codons
p_stop = regex.compile("\*")
p_x = regex.compile("X")

# Open fastq file(s) and search for patterns
for inFile in sys.argv[2:]:
    outFile = inFile + "-" + cellType + "-CDR3.csv"
    # extraFile = inFile + "-" + cellType + "-extra.txt"
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

    # try:
    #     fhExtra = open(extraFile, "w")
    # except:
    #     sys.exit("cannot write to file:" + extraFile)

    # Container to count stuff
    count_stuff = dict()

    for record in SeqIO.parse(fhIn, "fastq") :

        count_stuff["1. Total reads"] = count_stuff.get("1. Total reads",0) + 1

        # Translate sequence to protein
        translations = nucToPeptide(str(record.seq))

        # Find motif (V, J and everything in between)
        # TO DO: sometimes multiple patterns can match which results in a much longer CDR3. Needs a fix
        for i in range(len(translations)):
            # m = p.search(str(translations[i]))
            m = p.search(str(translations[i]))
            if m != None:   # print to file if there was a match

                # Extract CDR3 peptide sequence
                cdr3pep = m.group(1)
                aa_pos = list(m.span(1))

                # Remove the first 6, 5 or 4 aminoacids of the CDR3
                if cellType.startswith("IGH"):  # IGH
                    cdr3pep = cdr3pep[6:]
                    aa_pos[0] = aa_pos[0]+6
                elif cellType.startswith("IG"): # IGL, IGK
                    cdr3pep = cdr3pep[5:]
                    aa_pos[0] = aa_pos[0]+5
                else:                           # TRB, TRA
                    cdr3pep = cdr3pep[4:]
                    aa_pos[0] = aa_pos[0]+4

                # Remove the last 2 aminoacids of the CDR3
                cdr3pep = cdr3pep[:-2]
                aa_pos[1] = aa_pos[1] - 2

                # Check if there is an extra V motif, do this only with exact match of the V motif
                if p_v.search(cdr3pep) != None:
                    print("WARNING: Extra V motif in", record.id, cdr3pep)

                # Extract CDR3 nucleotide sequence
                nt_start = aa_pos[0] * 3
                nt_end = aa_pos[1] * 3
                if i < 3:                   # retrieve nucleotide reading frame +
                    tmp_seq = record.seq[i:]
                else:                       # retrieve nucleotide reading frame of reverse complement
                    tmp_seq = comrev(record.seq)[i-3:]
                cdr3nuc = tmp_seq[nt_start:nt_end]

                # Convert fastq quality scores to phred scores
                quality_scores = record.letter_annotations["phred_quality"]
                if i<3:
                    tmp_qual = quality_scores[i:]
                else:
                    tmp_qual = list(reversed(quality_scores))[i-3:]
                cdr3_quality_scores = tmp_qual[nt_start:nt_end]

                # Basic stats CDR3 quality
                cdr3_qual_min = min(cdr3_quality_scores)
                cdr3_qual_max = max(cdr3_quality_scores)
                cdr3_qual_avg = round(sum(cdr3_quality_scores)/float(len(cdr3_quality_scores)), 1)

                # Convert the quality scores to a string
                quality_scores = str(quality_scores)
                quality_scores = quality_scores.replace(",","").replace("[","").replace("]","")
                cdr3_quality_scores = str(cdr3_quality_scores)
                cdr3_quality_scores = cdr3_quality_scores.replace(",","").replace("[","").replace("]","")

                # Check for stop codons and untranslated codons
                if p_stop.search(cdr3pep) != None:    # stop codon in cdr3peptide
                    count_stuff["2. Discarded reads stop codon in CDR3"] = count_stuff.get("2. Discarded reads stop codon in CDR3",0) + 1
                    print("\t".join([record.id, str(i), str(record.seq), str(translations[i]), quality_scores]), file=fhStop)
                elif p_x.search(cdr3pep) != None:     # check for uncalled bases
                    count_stuff["3. Discarded reads uncalled bases in CDR3"] = count_stuff.get("3. Discarded reads uncalled bases in CDR3",0) + 1
                    print("\t".join([record.id, str(i), str(record.seq), str(translations[i]), quality_scores]), file=fhUncalled)
                else:                                 # Correct CDR3 peptide found without uncalled bases or stop codons
                    print("\t".join([record.id, str(i), str(cdr3pep), str(cdr3nuc), str(cdr3_qual_min), str(cdr3_qual_max), str(cdr3_qual_avg), cdr3_quality_scores]), file=fhOut)
                    print("\t".join([record.id, str(i), str(record.seq), str(translations[i]), quality_scores]), file=fhRaw)

                    # # Search for motif in the protein translation
                    # for m_extra in p_extra.finditer(str(translations[i])):
                    #     print(record.id, str(i), m_extra.group(0), m_extra.span(), file=fhExtra)

                    count_stuff["4. Reads with CDR3"] = count_stuff.get("4. Reads with CDR3",0) + 1
                    break
            else:
                print("\t".join([record.id, str(i), str(record.seq), str(translations[i])]), file=fhNoCdr3)

    # Make report
    print("Motifs:", motif, file=fhRep)

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
    # fhExtra.close()
