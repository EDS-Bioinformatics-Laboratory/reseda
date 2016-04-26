library(plyr)
data<-read.csv("/home/barbera/TMP/run06-clones_subs.csv", header=TRUE, sep="\t")
data$grp<-">1"
data$grp[which(data$read_perc<=1.0)] <- "0.9-1.0"
data$grp[which(data$read_perc<=0.9)] <- "0.8-0.9"
data$grp[which(data$read_perc<=0.8)] <- "0.7-0.8"
data$grp[which(data$read_perc<=0.7)] <- "0.6-0.7"
data$grp[which(data$read_perc<=0.6)] <- "0.5-0.6"
data$grp[which(data$read_perc<=0.5)] <- "0.4-0.5"
data$grp[which(data$read_perc<=0.4)] <- "0.3-0.4"
data$grp[which(data$read_perc<=0.3)] <- "0.2-0.3"
data$grp[which(data$read_perc<=0.2)] <- "0.1-0.2"
data$grp[which(data$read_perc<=0.1)] <- "0-0.1"

pivot<-ddply(data, .(grp), summarize, freq=length(grp))
total<-sum(pivot$freq)
pivot$perc<-100* pivot$freq / total
barplot(pivot$freq, names=pivot$grp, log="y", las=2, xlab="clonal size", ylab="frequency")
