from __future__ import print_function
from PairwiseAlignment import pairAlign # part of the lineage-tree repository
from sequences import readFasta


def alignAllVsAll(fileOut, sequences):
    '''
    Description: align all sequences against eachother
    In: fileOut, sequences[acc]=sequence
    Out: -, fileOut is written to disk
    '''

    fhOut = open(fileOut, "w")
    print("acc1 acc2 distance.all distance.replacements", file=fhOut)
    seqs = sequences.keys()   # [0:10]
    for i in range(0, len(seqs) - 1):
        idA = seqs[i]
        seqA = sequences[idA]

        for j in range(i + 1, len(seqs)):
            idB = seqs[j]
            seqB = sequences[idB]
            # print(idA, idB, seqA, seqB)
            print(idA, idB)
            aln = pairAlign(seqA, seqB)
            print(idA.split("|")[1], idB.split("|")[1], aln.distance, len(aln.pos_r), file=fhOut)
    fhOut.close()


def compareGenes(f):

    genes = readFasta(f)

    fileOut = f.replace(".fasta", ".distances.txt")
    alignAllVsAll(fileOut, genes)


if __name__ == '__main__':
    fastafiles = ["reference/TRBJ_human.fasta", "reference/TRBV_human.fasta", "reference/IGHJ_human.fasta", "reference/IGHV_human.fasta", "reference/TRBJ_mouse.fasta", "reference/TRBV_mouse.fasta", "reference/IGHJ_mouse.fasta", "reference/IGHV_mouse.fasta", "reference/IGKJ_human.fasta", "reference/IGKJ_mouse.fasta", "reference/IGKV_human.fasta", "reference/IGKV_mouse.fasta", "reference/IGLJ_human.fasta", "reference/IGLJ_mouse.fasta", "reference/IGLV_human.fasta", "reference/IGLV_mouse.fasta", "reference/TRAJ_human.fasta", "reference/TRAJ_mouse.fasta", "reference/TRAV_human.fasta", "reference/TRAV_mouse.fasta"]

    for f in fastafiles:
        print("Processing", f)
        compareGenes(f)

    print("FINISHED")
