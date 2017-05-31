V.file = "SP-Bsort14_S164_L001.assembled-CGTACGTA-IGHV_human-e-clean.sam.mut.txt"
J.file = "SP-Bsort14_S164_L001.assembled-CGTACGTA-IGHJ_human-e-clean.sam.mut.txt"
sample = gsub("(^.+)_L001.*", "\\1", V.file)
d.V = read.csv(V.file, header=T, sep=" ", stringsAsFactors = F)
d.J = read.csv(J.file, header=T, sep=" ", stringsAsFactors = F)

pdf(paste(sample, "-mutations.pdf", sep=""), width = 10, height = 10)
par(mfrow=c(2,2))
hist(d.V$mut.count, breaks=max(d.V$mut.count), col = "blue", main=paste("Mutations in V (count), n =",length(d.V[,1])))
hist(d.J$mut.count, breaks=max(d.J$mut.count), col = "green", main=paste("Mutations in J (count), n =",length(d.J[,1])))
hist(d.V$mut.perc, breaks=max(d.V$mut.perc), col = "blue", main=paste("Mutations/Length in V (%), n =",length(d.V[,1])))
hist(d.J$mut.perc, breaks=max(d.J$mut.perc), col = "green", main=paste("Mutations/Length in J (%), n =",length(d.J[,1])))
mtext(sample, side=3, outer=TRUE, line=-3)
dev.off()
