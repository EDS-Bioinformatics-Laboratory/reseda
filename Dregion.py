
# coding: utf-8

# In[1]:


import argparse
import pandas as pd
from Bio import SeqIO
from alignment.sequence import Sequence
from alignment.vocabulary import Vocabulary
from alignment.sequencealigner import SimpleScoring, LocalSequenceAligner
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


def pairLocalAlign(a, b, score_match, score_mismatch, score_gap):
        # Create sequences to be aligned.
    a = Sequence(a)
    b = Sequence(b)

    # Create a vocabulary and encode the sequences.
    v = Vocabulary()
    aEncoded = v.encodeSequence(a)
    bEncoded = v.encodeSequence(b)

    # Create a scoring and align the sequences using local aligner.
    scoring = SimpleScoring(score_match, score_mismatch)  # match, mismatch
    aligner = LocalSequenceAligner(scoring, score_gap)  # gap penalty, was -1
    score, encodeds = aligner.align(aEncoded, bEncoded, backtrace=True)

    # Iterate over optimal alignments and print them.
    for encoded in encodeds:
        alignment = v.decodeSequenceAlignment(encoded)

        # Get info about positions of (mis)matches and indels
        pos_m = list()
        pos_r = list()
        pos_i = list()
        pos_d = list()
        edit_transcript = list()
        for i in range(len(alignment.first)):
            if alignment.first[i] == alignment.second[i]:
                pos_m.append(i)
                edit_transcript.append("M")
            elif alignment.first[i] == "-":
                pos_i.append(i)
                edit_transcript.append("I")
            elif alignment.second[i] == "-":
                pos_d.append(i)
                edit_transcript.append("D")
            else:
                pos_r.append(i)
                edit_transcript.append("R")

        alignment.pos_m = pos_m
        alignment.pos_r = pos_r
        alignment.pos_i = pos_i
        alignment.pos_d = pos_d
        alignment.edit_transcript = edit_transcript

        # Edit distance
        alignment.distance = len(pos_r) + len(pos_d) + len(pos_i)

        return(alignment)


def getDeletedParts(refSeq, second):
    if "-" in second:
        alignparts = second.split("-")
        left = refSeq.split(alignparts[0])[0]
        right = refSeq.split(alignparts[-1])[-1]
        return(left, right)
    else:
        # MET HET ONDERSTAANDE: TOO MANY VALUES TO UNPACK: HIEROP TESTEN
        results = refSeq.split(second)
        if len(results) == 2:
            left, right = results
        else:
            left, right = ["unknown", "unknown"]
        return(left, right)


def createHistDens(pdf, df, column):
    # Density Plot and Histogram
    sns.distplot(df[column], hist=True, kde=True,
                 bins=int(180/5), color = 'darkblue',
                 hist_kws={'edgecolor':'black'},
                 kde_kws={'linewidth': 4})
    plt.title(column)
    pdf.savefig()
    plt.close()


# Input

parser = argparse.ArgumentParser(description='Detect the D region')
parser.add_argument('-cell', '--celltype', default='IGH', type=str, help='Cell type: (IGH|TRB|IGK|IGL|TRA) (default: %(default)s)')
parser.add_argument('-org', '--organism', default='human', type=str, help='Organism: (human|mouse) (default: %(default)s)')
parser.add_argument('-cdr3col', '--cdr3col', default='autodetect', type=str, help='Name of column that contains the CDR3 nucleotide sequence (default: %(default)s)')
parser.add_argument('-match', '--match', default=2, type=int, help='Alignment score for a match (default: %(default)s)')
parser.add_argument('-mismatch', '--mismatch', default=-1, type=int, help='Alignment score for a mismatch (default: %(default)s)')
parser.add_argument('-gap', '--gap', default=-5, type=int, help='Alignment penalty for a gap (default: %(default)s)')
parser.add_argument("cdr3_file", type=str, help='Path to cdr3 file. E.g. vjcdr3-clones-mut-PROJECT-IGH_HUMAN.csv')
args = parser.parse_args()

celltype = args.celltype
organism = args.organism
cdr3nuc_column = args.cdr3col
score_match = args.match
score_mismatch = args.mismatch
score_gap = args.gap

cdr3_file = args.cdr3_file

# Use this reference gene file
reference = celltype + "D" + "_" + organism + ".fasta"

# Output files
outfile = celltype + "D" + "_" + organism + "_" + str(score_match) + "_" + str(score_mismatch) + "_" + str(score_gap) + "_" + cdr3_file
outpdf = outfile.replace(".csv", ".pdf")
outsummary = outfile.replace(".csv", "-summary.csv")

# ## Read tabel containing a column with the CDR3 nucleotide
df = pd.read_csv(cdr3_file, sep="\t")

# Autodetect the cdr3nuc column or use user input
if cdr3nuc_column == "autodetect":
    if "cdr3nuc" in df.columns:
        cdr3nuc_column = "cdr3nuc"
    elif "cdr3nuc.mode" in df.columns:
        cdr3nuc_column = "cdr3nuc.mode"
    elif "cdr3nuc.min_mode" in df.columns:
        cdr3nuc_column = "cdr3nuc.min_mode"
    else:
        print("ERROR: couldn't automagically detect the column containing the CDR3 nucleotide sequence")
        exit()
    print("AUTO-DETECT CDR3 NUC COLUMN", cdr3nuc_column)
else:
    # use whatever the user specifies
    pass

# Get unique list of cdr3 nucleotide sequences to shorten time
cdr3nuc = df[cdr3nuc_column].unique().tolist()

# ## Read the reference sequences of the D region
try:
    fh = open(reference)
except:
    sys.exit("cannot open file: " + reference)

# Put sequences in a dictionary
sequences = dict()
for record in SeqIO.parse(fh, "fasta"):
    myId = str(record.id).split("|")[1]  # in case of IMGT reference sequences
    sequences[myId] = record.seq

# Perform pairwise sequence alignment between all CDR3's and all reference D genes
cdr3nuc_list = list()
idD_list = list()
seqD_list = list()
alignFirst = list()
alignSecond = list()
edit_transcript = list()
distance = list()
gapCount = list()
identicalCount = list()
score = list()
similarCount = list()

for seqA in cdr3nuc:
    if pd.isna(seqA): # Skip if value is empty
        continue
    seqA = seqA.upper()
    for refID, seqB in sequences.items():
        seqB = seqB.upper()
        alignment = pairLocalAlign(seqA, seqB, score_match, score_mismatch, score_gap)

        cdr3nuc_list.append(seqA)
        idD_list.append(refID)
        seqD_list.append("".join(seqB))
        alignFirst.append("".join(alignment.first))
        alignSecond.append("".join(alignment.second))
        edit_transcript.append("".join(alignment.edit_transcript))
        distance.append(alignment.distance)
        gapCount.append(alignment.gapCount)
        identicalCount.append(alignment.identicalCount)
        score.append(alignment.score)
        similarCount.append(alignment.similarCount)

df_align = pd.DataFrame({'cdr3nuc': cdr3nuc_list, 'refID': idD_list, 'refSeq': seqD_list,
                         'first': alignFirst, 'second': alignSecond, 'edit_transcript': edit_transcript,
                         'distance': distance, 'gapCount': gapCount, 'identicalCount': identicalCount,
                         'score': score, 'similarCount': similarCount})

df_align = df_align.sort_values(by=['cdr3nuc', 'score'], ascending=False)

# ## Get highest scoring matches per CDR3
df_high_score = df_align.groupby('cdr3nuc').agg({'score': max})
df_high_score = df_high_score.reset_index()

df_best_align = pd.merge(df_align, df_high_score, on=['cdr3nuc', 'score'])

# ## What kind of measures can we use to check the performance?
#
# * nr of matches and mismatches (match should be high, mismatches are expected but should be lower
# * nr of indels (should be none or low)
# * nr of nucleotides that are deleted from the ends
# * length of the alignment
#
# Do this on BCRh and on TCRb. The latter doesn't have SHM, so there shouldn't be mismatches, except for some sequencing errors

f = lambda x: x.count("M")
df_best_align["countM"] = [x for x in map(f, df_best_align["edit_transcript"])]
f = lambda x: x.count("R")
df_best_align["countR"] = [x for x in map(f, df_best_align["edit_transcript"])]
f = lambda x: x.count("I")
df_best_align["countI"] = [x for x in map(f, df_best_align["edit_transcript"])]
f = lambda x: x.count("D")
df_best_align["countD"] = [x for x in map(f, df_best_align["edit_transcript"])]
df_best_align["countIndels"] = df_best_align["countI"] + df_best_align["countD"]

f = lambda x: len(x)
df_best_align["alignLength"] = [x for x in map(f, df_best_align["edit_transcript"])]

# Determine the parts of the reference D gene that was deleted
results = [x for x in map(getDeletedParts, df_best_align["refSeq"], df_best_align["second"])]
returnLeft = lambda x: x[0]
returnRight = lambda x: x[1]
df_best_align["deletedLeft"] = [x for x in map(returnLeft, results)]
df_best_align["deletedRight"] = [x for x in map(returnRight, results)]

# Lengths of the deleted nucleotides
delLength = lambda x: len(x)
df_best_align["countDeletedLeft"] = [x for x in map(delLength, df_best_align["deletedLeft"])]
df_best_align["countDeletedRight"] = [x for x in map(delLength, df_best_align["deletedRight"])]

# Make histograms and density plots
cols = ['distance', 'gapCount', 'identicalCount', 'score', 'similarCount',
       'countM', 'countR', 'countI', 'countD', 'countIndels', 'alignLength',
       'countDeletedLeft', 'countDeletedRight']

#with PdfPages(outpdf) as pdf:
#    for col in cols:
#        createHistDens(pdf, df_best_align, col)
#print("Wrote", outpdf, "to disk")

# Write the results to file
df_best_align.to_csv(outfile)
print("Wrote", outfile, "to disk")

# Calculate mean, median, standard deviation, outlier? and store in table
# Write this table to disk
means = list()
medians = list()
stddevs = list()
for col in cols:
    means.append(df_best_align[col].mean())
    medians.append(df_best_align[col].median())
    stddevs.append(df_best_align[col].std())
df_summary = pd.DataFrame({'align_info': cols, 'mean': means, 'median': medians, 'stddev': stddevs})
df_summary["organism"] = organism
df_summary["celltype"] = celltype
df_summary["score_match"] = score_match
df_summary["score_mismatch"] = score_mismatch
df_summary["score_gap"] = score_gap
df_summary["cdr3file"] = cdr3_file
df_summary["cdr3nuc_column"] = cdr3nuc_column

df_summary.to_csv(outsummary)
print("Wrote", outsummary, "to disk")
