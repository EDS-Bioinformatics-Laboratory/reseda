library(plyr)

f1="/mnt/immunogenomics/RUNS/run06-20160306-miseq/results-tbcell/run-clones_subs-TRB_HUMAN-after-reassignment.csv"
f2="/mnt/immunogenomics/RUNS/run06-20160306-miseq/results-tbcell-nov2016/run-clones_subs-TRB_HUMAN-after-reassignment.csv"

clones1=read.csv(f1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
clones2=read.csv(f2, header=TRUE, sep="\t", stringsAsFactors=FALSE)
both=merge(clones1,clones2, by=c("Sample","V_sub","J_sub","cdr3pep"), all=TRUE)
both$freq.x[which(is.na(both$freq.x))] = 0.1
both$freq.y[which(is.na(both$freq.y))] = 0.1

plot(both$freq.x, both$freq.y, log="xy", main="Clones (VJCDR3)", xlab="Without cdr3_min_qual filter", ylab="With cdr3_min_qual filter")

high.clones1 = both[which(both$freq.y==0.1),]
high.clones1 = high.clones1[order(high.clones1$freq.x, decreasing = TRUE),]
high.clones2 = both[which(both$freq.x==0.1),]
high.clones2 = high.clones2[order(high.clones2$freq.y, decreasing = TRUE),]

head(high.clones1)
head(high.clones2)
