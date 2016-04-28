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

pivot<-ddply(data, .(Sample,grp), summarize, freq=length(grp))
pivot<-pivot[order(pivot$grp,pivot$Sample),]
barplot(pivot$freq, names=pivot$grp, beside=TRUE, log="y", las=2, xlab="clonal size", ylab="average frequency")


avg_pivot<-ddply(pivot, .(grp), summarize, avg=mean(freq), sd=sd(freq), n=length(grp))
avg_pivot$se <- avg_pivot$sd / sqrt(avg_pivot$n)

barCenters<-barplot(avg_pivot$avg, names=avg_pivot$grp, beside=TRUE, log="y", las=2, xlab="clonal size", ylab="frequency")

# Add error bars. See: http://datascienceplus.com/building-barplots-with-error-bars/
arrows(barCenters, avg_pivot$avg - avg_pivot$se * 2, barCenters,
       avg_pivot$avg + avg_pivot$se * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

# Output table as latex
library(xtable)
xtable(avg_pivot)
