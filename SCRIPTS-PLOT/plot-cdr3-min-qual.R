file.name<-"SP-CB28_S61_L001.assembled-CATGCATG-IGH_HUMAN-all_info.csv.min_qual.csv"
data<-read.csv(file.name, sep=" ")
names(data)<-c("min.qual","freq")

# Calculate percentage
total<-sum(data$freq)
data$perc<-data$freq/total

# Plot percentage of reads with minimum quality value in CDR3
png(file = paste(file.name,".png",sep=""), width=1200, height=800)
barplot(data$perc, names.arg = data$min.qual, main=file.name, xlab="Minimum base quality in CDR3", ylab="Percentage of reads")
dev.off()
