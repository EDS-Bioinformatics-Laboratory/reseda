plotMutations<-function(directory, V.file, J.file){
  sample = gsub("(^.+)_L001.*", "\\1", V.file)
  d.V = read.csv(paste(directory,V.file,sep=""), header=T, sep=" ", stringsAsFactors = F)
  d.J = read.csv(paste(directory,J.file,sep=""), header=T, sep=" ", stringsAsFactors = F)
  
  # Remove sequences without mutations
  #d.V = d.V[which(d.V$mut.count>0),]
  #d.J = d.J[which(d.J$mut.count>0),]
  
  # Plot
  pdf(paste(sample, "-mutations.pdf", sep=""), width = 10, height = 10)
  par(mfrow=c(2,2))
  hist(d.V$mut.count, breaks=max(d.V$mut.count), col = "blue", main=paste("Mutations in V (count), n =",length(d.V[,1])))
  hist(d.J$mut.count, breaks=max(d.J$mut.count), col = "green", main=paste("Mutations in J (count), n =",length(d.J[,1])))
  hist(d.V$mut.frac, breaks=20, col = "blue", main=paste("Mutations/Length in V, n =",length(d.V[,1])))
  hist(d.J$mut.frac, breaks=20, col = "green", main=paste("Mutations/Length in J, n =",length(d.J[,1])))
  mtext(sample, side=3, outer=TRUE, line=-3)
  dev.off()
  cat("Wrote", paste(sample, "-mutations.pdf", sep=""), "to disk")
}

# Read data
directory = "/mnt/immunogenomics/RUNS/run07-20160401-miseq/results-tbcell/raw/"
V.file = "SP-Bsort06_S157_L001.assembled-ACGTACGT-IGHV_human-e-clean.sam.mut.txt" # Naive
J.file = "SP-Bsort06_S157_L001.assembled-ACGTACGT-IGHJ_human-e-clean.sam.mut.txt"
plotMutations(directory, V.file, J.file)

V.file = "SP-Bsort07_S158_L001.assembled-ACTGACTG-IGHV_human-e-clean.sam.mut.txt" # Memory
J.file = "SP-Bsort07_S158_L001.assembled-ACTGACTG-IGHJ_human-e-clean.sam.mut.txt"
plotMutations(directory, V.file, J.file)

V.file = "SP-Bsort09_S160_L001.assembled-ATCGATCG-IGHV_human-e-clean.sam.mut.txt" # Paxgene
J.file = "SP-Bsort09_S160_L001.assembled-ATCGATCG-IGHJ_human-e-clean.sam.mut.txt"
plotMutations(directory, V.file, J.file)
