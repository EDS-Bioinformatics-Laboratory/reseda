library("plyr")

# Read CDR3 lengths
cdr3.length.all=read.csv("run06-tcrb-cdr3-length-all-but-reassigned.txt", header = FALSE)
cdr3.length.reassigned=read.csv("run06-tcrb-cdr3-length-reassigned.txt", header = FALSE)
colnames(cdr3.length.all)=c("CDR3length")
colnames(cdr3.length.reassigned)=c("CDR3length")

# Histogram of the CDR3 length
hist.cdr3.length.all=ddply(cdr3.length.all, .(CDR3length), summarize, freq=length(CDR3length))
hist.cdr3.length.reassigned=ddply(cdr3.length.reassigned, .(CDR3length), summarize, freq=length(CDR3length))
plot(hist.cdr3.length.all, type="o", col="blue")
lines(hist.cdr3.length.reassigned, type="o", pch=22, lty=2, col="red")

# Combine in one table and calculate fraction of re-assigned
hist.combined=merge(hist.cdr3.length.all, hist.cdr3.length.reassigned, by="CDR3length", all=TRUE)
colnames(hist.combined) = c("CDR3length","all","reassigned")
#hist.combined$reassigned[which(is.na(hist.combined$reassigned))]=0
hist.combined$frac.reassigned=hist.combined$reassigned/(hist.combined$reassigned+hist.combined$all)

# Discard non info after length 31
hist.trunc=hist.combined[1:28,]

# Plot re-assigned vs CDR3 length
plot(hist.trunc$CDR3length,hist.trunc$frac.reassigned, type="o", col="purple", xlab="CDR3 length", ylab="Fraction re-assigned")
