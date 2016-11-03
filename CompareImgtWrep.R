library(plyr)

setwd("/home/barbera/Data/tbcell/RepSeq2016")

imgt=read.csv("Compare-WREP-IMGT/1_Summary.txt.parsed.csv", header=TRUE, sep=" ", stringsAsFactors=FALSE)
wrep=read.csv("Compare-WREP-IMGT/S074-057_S102_L001.assembled-CATGCATG-TRB_HUMAN-all_info.csv.rr.all_info.csv.parsed.csv", header=TRUE, sep=" ", stringsAsFactors=FALSE)
both=merge(imgt,wrep, by="acc")

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

both.unequal.cdr3 = both[which(both$cdr3.x != both$cdr3.y),]
both.unequal.v = both[which(both$v.x != both$v.y),]
both.unequal.j = both[which(both$j.x != both$j.y),]
both.unequal.vjcdr3 = both[which(both$cdr3.x != both$cdr3.y | both$v.x != both$v.y | both$j.x != both$j.y),]

# Clones

imgt.clones.cdr3 = ddply(imgt, .(cdr3), summarize, freq=length(acc))
wrep.clones.cdr3 = ddply(wrep, .(cdr3), summarize, freq=length(acc))
imgt.clones.vjcdr3 = ddply(imgt, .(v,j,cdr3), summarize, freq=length(acc))
wrep.clones.vjcdr3 = ddply(wrep, .(v,j,cdr3), summarize, freq=length(acc))

clones.cdr3 = merge(imgt.clones.cdr3, wrep.clones.cdr3, by="cdr3", all=TRUE)
clones.cdr3$freq.x[which(is.na(clones.cdr3$freq.x))] = 0.1
clones.cdr3$freq.y[which(is.na(clones.cdr3$freq.y))] = 0.1
plot(clones.cdr3$freq.x, clones.cdr3$freq.y, log="xy", main="Clones (CDR3)", xlab="IMGT", ylab="WREP")

clones.vjcdr3 = merge(imgt.clones.vjcdr3, wrep.clones.vjcdr3, by=c("v","j","cdr3"), all=TRUE)
clones.vjcdr3$freq.x[which(is.na(clones.vjcdr3$freq.x))] = 0.1
clones.vjcdr3$freq.y[which(is.na(clones.vjcdr3$freq.y))] = 0.1
plot(clones.vjcdr3$freq.x, clones.vjcdr3$freq.y, log="xy", main="Clones (VJCDR3)", xlab="IMGT", ylab="WREP")

high.imgt.clones.vjcdr3 = clones.vjcdr3[which(clones.vjcdr3$freq.y == 0.1),]
high.imgt.clones.vjcdr3 = high.imgt.clones.vjcdr3[order(high.imgt.clones.vjcdr3$freq.x, decreasing=TRUE),]
high.wrep.clones.vjcdr3 = clones.vjcdr3[which(clones.vjcdr3$freq.x == 0.1),]
high.wrep.clones.vjcdr3 = high.wrep.clones.vjcdr3[order(high.wrep.clones.vjcdr3$freq.y, decreasing=TRUE),]



head(clones.cdr3)
plot(clones.cdr3$freq.x, clones.cdr3$freq.y)
plot(clones.cdr3$freq.x, clones.cdr3$freq.y, log=TRUE)
plot(clones.cdr3$freq.x, clones.cdr3$freq.y, log="xy")
plot(clones.cdr3$freq.x, clones.cdr3$freq.y, log="xy", main="Clones (CDR3)", xlab="IMGT", ylab="WREP")
head(clones.cdr3)
is.na(clones.cdr3$freq.x)
which(clones.freq.x==NA)
which(clones$freq.x==NA)
which(clones.cdr3$freq.x==NA)
which(is.na(clones.cdr3$freq.x))
plot(clones.cdr3$freq.x, clones.cdr3$freq.y, log="xy", main="Clones (CDR3)", xlab="IMGT", ylab="WREP")
clones.cdr3$freq.x[which(is.na(clones.cdr3$freq.x))] = 0
head(clones.cdr3)
clones.cdr3$freq.y[which(is.na(clones.cdr3$freq.y))] = 0
head(clones.cdr3)
plot(clones.cdr3$freq.x, clones.cdr3$freq.y, log="xy", main="Clones (CDR3)", xlab="IMGT", ylab="WREP")
clones.cdr3 = merge(imgt.clones.cdr3, wrep.clones.cdr3, by="cdr3", all=TRUE)
clones.cdr3$freq.x[which(is.na(clones.cdr3$freq.x))] = 0.1
clones.cdr3$freq.y[which(is.na(clones.cdr3$freq.y))] = 0.1
plot(clones.cdr3$freq.x, clones.cdr3$freq.y, log="xy", main="Clones (CDR3)", xlab="IMGT", ylab="WREP")
head(imgt)
imgt.clones.vjcdr3 = ddply(imgt, .(v,j,cdr3), summarize, freq=length(acc))
wrep.clones.vjcdr3 = ddply(wrep, .(v,j,cdr3), summarize, freq=length(acc))
clones.vjcdr3 = merge(imgt.clones.vjcdr3, wrep.clones.vjcdr3, by=c("v","j","cdr3"), all=TRUE)
head(clones.vjcdr3)
clones.vjcdr3$freq.x[which(is.na(clones.vjcdr3$freq.x))] = 0.1
clones.vjcdr3$freq.y[which(is.na(clones.vjcdr3$freq.y))] = 0.1
plot(clones.vjcdr3$freq.x, clones.vjcdr3$freq.y, log="xy", main="Clones (VJCDR3)", xlab="IMGT", ylab="WREP")
which(clones.cdr3$freq.x==0.1)
head(clones.cdr3[which(clones.cdr3$freq.x==0.1),])
head(clones.cdr3[which(clones.cdr3$freq.x==0.1),], n = 50)
head(clones.cdr3[which(clones.cdr3$freq.y==0.1),], n = 50)
View(both.unequal.cdr3)
View(both.unequal.v)
which(both.unequal.v$v.y=="TRBV12-3+TRVB12-4")
which(both.unequal.v$v.y=="TRBV12-3+TRBV12-4")
length(which(both.unequal.v$v.y=="TRBV12-3+TRBV12-4"))
length(which(both.unequal.v$v.x=="TRBV12-3"))
length(which(both.unequal.v$v.x=="TRBV12-4"))
length(which(both.unequal.v$v.y=="TRBV6-5+TRBV6-6"))
cat("CDR3 OVERLAP")
cat("Equal CDR3:", cdr3.equal, 100*cdr3.equal/total.acc)
cat("Unequal CDR3:", cdr3.unequal, 100*cdr3.unequal/total.acc)
cat("V GENE OVERLAP")
cat("Equal V gene:", v.equal, 100*v.equal/total.acc)
cat("Unequal V gene:", v.unequal, 100*v.unequal/total.acc)
cat("J GENE OVERLAP")
cat("Equal J gene:", j.equal, 100*j.equal/total.acc)
cat("Unequal J gene:", j.unequal, 100*j.unequal/total.acc)
both.unequal.vjcdr3 = both[which(both$cdr3.x != both$cdr3.y || both$v.x != both$v.y || both$j.x != both$j.y)]
both.unequal.vjcdr3 = both[which(both$cdr3.x != both$cdr3.y || both$v.x != both$v.y || both$j.x != both$j.y),]
both.unequal.vjcdr3 = both[which(both$cdr3.x != both$cdr3.y | both$v.x != both$v.y | both$j.x != both$j.y),]
View(both.unequal.vjcdr3)
high.imgt.clones.vjcdr3 = clones.vjcdr3[which(clones.vjcdr3$freq.y == 0.1),]
order(clones.vjcdr3$freq.x, desc=TRUE)
help(order)
order(clones.vjcdr3$freq.x, descreasing=TRUE)
order(clones.vjcdr3$freq.x, decreasing=TRUE)
high.imgt.clones.vjcdr3 = high.imgt.clones.vjcdr3[order(clones.vjcdr3$freq.x, decreasing=TRUE),]
head(high.imgt.clones.vjcdr3)
high.imgt.clones.vjcdr3 = clones.vjcdr3[which(clones.vjcdr3$freq.y == 0.1),]
high.imgt.clones.vjcdr3 = high.imgt.clones.vjcdr3[order(high.imgt.clones.vjcdr3$freq.x, decreasing=TRUE),]
View(high.imgt.clones.vjcdr3)
high.wrep.clones.vjcdr3 = clones.vjcdr3[which(clones.vjcdr3$freq.x == 0.1),]
high.wrep.clones.vjcdr3 = high.wrep.clones.vjcdr3[order(high.wrep.clones.vjcdr3$freq.y, decreasing=TRUE),]
View(high.wrep.clones.vjcdr3)
View(both.unequal.cdr3)