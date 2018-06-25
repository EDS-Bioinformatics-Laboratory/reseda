from __future__ import print_function
import sys
import argparse
from Bio import SeqIO


def getGeneNames(ref):
    '''
    Description: open fasta file, make dictionary d[accession] = fasta_description
    In: fasta file name (str)
    Out: d[accession] = fasta_description (dict)
    '''
    try:
        fh = open(ref)
    except:
        sys.exit("Cannot open file: " + ref)

    d = dict()
    d["*"] = "*"
    for line in fh:
        if line.startswith('>'):
            line = line.rstrip()
            line = line.split(" ")
            acc = line[0]
            descr = " ".join(line[1:])
            d[acc] = descr
    fh.close()
    return(d)


def summaryGeneNames(refseq, samfile):
    '''
    Description: lookup gene name for accessions in sam file
    In: refseq (dict), sam file name (str)
    Out: file name with gene counts (str)
    '''
    try:
        fh = open(samfile)
    except:
        sys.exit("cannot open file: " + samfile)

    gene_summary = dict()  # to count gene names
    for line in fh:
        if line.startswith("@"): # skip header
            continue
        line = line.rstrip()
        line = line.split("\t")
        acc = line[2]
        gene = refseq.get(acc, "UNKNOWN")
        gene_summary[gene] = gene_summary.get(gene, 0) + 1
    fh.close()

    outfile = sam_files + "-genes.txt"
    try:
        fhOut = open(outfile, "w")
    except:
        sys.exit("cannot write file: " + outfile)

    for gene in sorted(gene_summary, key=gene_summary.get, reverse=True):
        print(gene, gene_summary[gene], file=fhOut)
    fhOut.close()
    return(outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Open SAM file and print the RefSeq gene names')
    parser.add_argument('-ref', '--reference', default='reference/GRCh38_latest_rna.fna', type=str, help='Fasta file of the RefSeq database (default: %(default)s)')
    parser.add_argument("sam_files", type=str, nargs='+', help='Path(s) to SAM file(s)')
    args = parser.parse_args()

    # make dictionary of acc -> fasta_description
    refseq = getGeneNames(args.reference)

    # Open each sam file, lookup gene name and report
    for samfile in args.sam_files:
        outfile = summaryGeneNames(refseq, samfile)
        print("Wrote", outfile, "to disk")

    print("FINISHED")
