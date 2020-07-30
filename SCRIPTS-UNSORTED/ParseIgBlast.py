import sys
import argparse

def parseEntry(entry):
    lines = entry.split("\n")
    for i, line in enumerate(lines):
        if line.startswith("# Query: "):
            query = line.replace("# Query: ", "")
            print(query)
        if line.startswith("# V-(D)-J rearrangement summary"):
            rearrangementinfo = line[i+1]
        if line.startswith("# V-(D)-J junction details"):
            junctioninfo = line[i+1]
        if line.startswith("# Sub-region sequence details"):
            cdr3info = line[i+1]
        if line.startswith("# Alignment summary"):
            pass   # read all lines up to the empty line (it contains the FWR and CDR locations)

def parseIgBlast(blast_file, entry_delim="# IGBLASTN"):
    fh = open(blast_file)
    entry = ""
    for line in fh:
        if line.startswith(entry_delim) and entry != "":
            parseEntry(entry)
            entry = ""
        entry = entry + line
    parseEntry(entry) # do not forget to process the last entry


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse IgBlast file (run it like this: ./bin/igblastn -germline_db_V database/IGHV_human_igBlast.fasta -germline_db_D database/IGHD_human_igBlast.fasta -germline_db_J database/IGHJ_human_igBlast.fasta -query ../test-run32-B_S58_L001.assembled.fasta -auxiliary_data optional_file/human_gl.aux -show_translation -outfmt '7 std qseq sseq btop' > blast2.out)")
    #parser.add_argument('-i', '--input', default='blast.out', type=str, help='IgBlast output file (default: %(default)s)')
    parser.add_argument("blast_files", type=str, nargs='+', help='Path(s) to IgBlast output file(s)')
    args = parser.parse_args()

    for blast_file in args.blast_files:
        parseIgBlast(blast_file)
