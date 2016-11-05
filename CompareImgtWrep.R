library(plyr)

#setwd("/home/barbera/Data/tbcell/RepSeq2016/Compare-WREP-IMGT")
setwd("/home/narya/TMP/Compare-WREP-IMGT")

imgt=read.csv("1_Summary.txt.parsed.csv", header=TRUE, sep=" ", stringsAsFactors=FALSE)
wrep=read.csv("S074-057_S102_L001.assembled-CATGCATG-TRB_HUMAN-all_info.csv.rr.all_info.csv.parsed.csv", header=TRUE, sep=" ", stringsAsFactors=FALSE)
both=merge(imgt,wrep, by="acc")

# Check how much overlap there is based on CDR3, V gene or J gene
total.acc = length(both[,1])
cdr3.equal = length(which(both$cdr3.x == both$cdr3.y))
cdr3.unequal = length(which(both$cdr3.x != both$cdr3.y))
v.equal = length(which(both$v.x == both$v.y))
v.unequal = length(which(both$v.x != both$v.y))
j.equal = length(which(both$j.x == both$j.y))
j.unequal = length(which(both$j.x != both$j.y))

cat("CDR3 OVERLAP")
cat("Equal CDR3:", cdr3.equal, 100*cdr3.equal/total.acc)
cat("Unequal CDR3:", cdr3.unequal, 100*cdr3.unequal/total.acc)
cat("V GENE OVERLAP")
cat("Equal V gene:", v.equal, 100*v.equal/total.acc)
cat("Unequal V gene:", v.unequal, 100*v.unequal/total.acc)
cat("J GENE OVERLAP")
cat("Equal J gene:", j.equal, 100*j.equal/total.acc)
cat("Unequal J gene:", j.unequal, 100*j.unequal/total.acc)

# Where are the differences in CDR3, V, J or the combination
both.unequal.cdr3 = both[which(both$cdr3.x != both$cdr3.y),]
both.unequal.v = both[which(both$v.x != both$v.y),]
both.unequal.j = both[which(both$j.x != both$j.y),]
both.unequal.vjcdr3 = both[which(both$cdr3.x != both$cdr3.y | both$v.x != both$v.y | both$j.x != both$j.y),]

# On the clone level

imgt.clones.cdr3 = ddply(imgt, .(cdr3), summarize, freq=length(acc))
wrep.clones.cdr3 = ddply(wrep, .(cdr3), summarize, freq=length(acc))
imgt.clones.vjcdr3 = ddply(imgt, .(v,j,cdr3), summarize, freq=length(acc))
wrep.clones.vjcdr3 = ddply(wrep, .(v,j,cdr3), summarize, freq=length(acc))

# Combine data on CDR3
clones.cdr3 = merge(imgt.clones.cdr3, wrep.clones.cdr3, by="cdr3", all=TRUE)
clones.cdr3$freq.x[which(is.na(clones.cdr3$freq.x))] = 0.1
clones.cdr3$freq.y[which(is.na(clones.cdr3$freq.y))] = 0.1
plot(clones.cdr3$freq.x, clones.cdr3$freq.y, log="xy", main="Clones (CDR3)", xlab="IMGT", ylab="WREP")

# Combine data on VJCDR3
clones.vjcdr3 = merge(imgt.clones.vjcdr3, wrep.clones.vjcdr3, by=c("v","j","cdr3"), all=TRUE)
clones.vjcdr3$freq.x[which(is.na(clones.vjcdr3$freq.x))] = 0.1
clones.vjcdr3$freq.y[which(is.na(clones.vjcdr3$freq.y))] = 0.1
plot(clones.vjcdr3$freq.x, clones.vjcdr3$freq.y, log="xy", main="Clones (VJCDR3)", xlab="IMGT", ylab="WREP")

# Get clones that are in one dataset, sort on clone frequency
high.imgt.clones.vjcdr3 = clones.vjcdr3[which(clones.vjcdr3$freq.y == 0.1),]
high.imgt.clones.vjcdr3 = high.imgt.clones.vjcdr3[order(high.imgt.clones.vjcdr3$freq.x, decreasing=TRUE),]
high.wrep.clones.vjcdr3 = clones.vjcdr3[which(clones.vjcdr3$freq.x == 0.1),]
high.wrep.clones.vjcdr3 = high.wrep.clones.vjcdr3[order(high.wrep.clones.vjcdr3$freq.y, decreasing=TRUE),]

# Count overlap on the clone level (vjcdr3)
total.clones=length(clones.vjcdr3[,1])
clones.only.in.imgt=length(high.imgt.clones.vjcdr3[,1])
clones.only.in.wrep=length(high.wrep.clones.vjcdr3[,1])
cat("Total clones", total.clones)
cat("Clones only in IMGT", clones.only.in.imgt, 100*clones.only.in.imgt/total.clones)
cat("Clones only in WREP", clones.only.in.wrep, 100*clones.only.in.wrep/total.clones)

# Show and compare the top100 clones that are in one dataset
top100.imgt.clones.vjcdr3 = head(high.imgt.clones.vjcdr3, n=100)
top100.wrep.clones.vjcdr3 = head(high.wrep.clones.vjcdr3, n=100)
top100.compare = merge(top100.imgt.clones.vjcdr3, top100.wrep.clones.vjcdr3, by="cdr3", all=TRUE)
top100.compare = top100.compare[order(top100.compare$freq.x.x, decreasing = TRUE),]
write.csv(top100.compare, "top100.compare.csv")
