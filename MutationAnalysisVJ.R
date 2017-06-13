library(plyr)

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
  cat("Wrote", paste(sample, "-mutations.pdf", sep=""), "to disk\n")
  
  return(list(d.V, d.J, sample))
}

main<-function(directory, V.file, J.file, CDR3.file){
  # Read and combine data
  d = plotMutations(directory, V.file, J.file)
  d.V = d[[1]]
  d.J = d[[2]]
  sample = d[[3]]
  rm(d)  # clean up
  d.CDR3 = read.csv(CDR3.file, header=T, sep="\t", stringsAsFactors = F)
  d.CDR3$VJCDR3 = paste(d.CDR3$V_sub, d.CDR3$J_sub, d.CDR3$cdr3pep, sep = "-")
  d.combined = merge(merge(d.CDR3, d.V, by="acc"), d.J, by="acc")
  write.csv(d.combined, file=paste(sample, "-mutations-per-accession.csv", sep=""))
  cat("Wrote", paste(sample,"-mutations-per-accession.csv",sep=""), "to disk\n")
  
  # Summarize mutations per VJCDR3 clone
  clones = ddply(d.combined, .(VJCDR3), summarise, freq=length(acc), total.mut.count.V=sum(mut.count.x), avg.mut.frac.V=mean(mut.frac.x), total.mut.count.J=sum(mut.count.y), avg.mut.frac.J=mean(mut.frac.y))
  clones = clones[order(clones$freq, decreasing=T),]
  write.csv(clones,file=paste(sample,"-mutations-per-clone.csv",sep=""))
  cat("Wrote", paste(sample,"-mutations-per-clone.csv",sep=""), "to disk\n")
  return(clones)
}

# Read data
directory = "/mnt/immunogenomics/RUNS/run07-20160401-miseq/results-tbcell/raw/"
V.file = "SP-Bsort06_S157_L001.assembled-ACGTACGT-IGHV_human-e-clean.sam.mut.txt" # Naive
J.file = "SP-Bsort06_S157_L001.assembled-ACGTACGT-IGHJ_human-e-clean.sam.mut.txt"
CDR3.file = "/mnt/immunogenomics/RUNS/run07-20160401-miseq/results-tbcell/final/correct-mid/SP-Bsort06_S157_L001.assembled-ACGTACGT-IGH_HUMAN-all_info.csv.rr.all_info.csv"
clones=main(directory, V.file, J.file, CDR3.file)

V.file = "SP-Bsort07_S158_L001.assembled-ACTGACTG-IGHV_human-e-clean.sam.mut.txt" # Memory
J.file = "SP-Bsort07_S158_L001.assembled-ACTGACTG-IGHJ_human-e-clean.sam.mut.txt"
CDR3.file = "/mnt/immunogenomics/RUNS/run07-20160401-miseq/results-tbcell/final/correct-mid/SP-Bsort07_S158_L001.assembled-ACTGACTG-IGH_HUMAN-all_info.csv.rr.all_info.csv"
clones=main(directory, V.file, J.file, CDR3.file)

V.file = "SP-Bsort09_S160_L001.assembled-ATCGATCG-IGHV_human-e-clean.sam.mut.txt" # Paxgene
J.file = "SP-Bsort09_S160_L001.assembled-ATCGATCG-IGHJ_human-e-clean.sam.mut.txt"
CDR3.file = "/mnt/immunogenomics/RUNS/run07-20160401-miseq/results-tbcell/final/correct-mid/SP-Bsort09_S160_L001.assembled-ATCGATCG-IGH_HUMAN-all_info.csv.rr.all_info.csv"
clones=main(directory, V.file, J.file, CDR3.file)
