setwd("/home/barbera/git/lineage-tree/")

t10=read.csv("AB-RBF126-T_S92_L001.assembled-ATCGATCG-TRB_HUMAN-all_info.csv.rr.all_info.csv.report.10.txt", header=T, sep=" ", stringsAsFactors = F)
t20=read.csv("AB-RBF126-T_S92_L001.assembled-ATCGATCG-TRB_HUMAN-all_info.csv.rr.all_info.csv.report.20.txt", header=T, sep=" ", stringsAsFactors = F)
t30=read.csv("AB-RBF126-T_S92_L001.assembled-ATCGATCG-TRB_HUMAN-all_info.csv.rr.all_info.csv.report.30.txt", header=T, sep=" ", stringsAsFactors = F)
t40=read.csv("AB-RBF126-T_S92_L001.assembled-ATCGATCG-TRB_HUMAN-all_info.csv.rr.all_info.csv.report.40.txt", header=T, sep=" ", stringsAsFactors = F)

t10=read.csv("AB-RBF147-B_S143_L001.assembled-ACTGACTG-IGH_HUMAN-all_info.csv.rr.all_info.csv.report.10.txt", header=T, sep=" ", stringsAsFactors = F)
t20=read.csv("AB-RBF147-B_S143_L001.assembled-ACTGACTG-IGH_HUMAN-all_info.csv.rr.all_info.csv.report.20.txt", header=T, sep=" ", stringsAsFactors = F)
t30=read.csv("AB-RBF147-B_S143_L001.assembled-ACTGACTG-IGH_HUMAN-all_info.csv.rr.all_info.csv.report.30.txt", header=T, sep=" ", stringsAsFactors = F)
t40=read.csv("AB-RBF147-B_S143_L001.assembled-ACTGACTG-IGH_HUMAN-all_info.csv.rr.all_info.csv.report.40.txt", header=T, sep=" ", stringsAsFactors = F)

boxplot(t10$low_qual, t20$low_qual, t30$low_qual, t40$low_qual, names=c("10", "20", "30", "40"), xlab="Base quality threshold", ylab="Read count low quality CDR3")
boxplot(t10$high_qual, t20$high_qual, t30$high_qual, t40$low_qual, names=c("10", "20", "30", "40"), xlab="Base quality threshold", ylab="Read count high quality CDR3")

plot(t10$low_qual, t10$high_qual, main=paste("n=",length(t10$CDR3)))
plot(t20$low_qual, t20$high_qual, main=paste("n=",length(t20$CDR3)))
plot(t30$low_qual, t30$high_qual, main=paste("n=",length(t30$CDR3)))
plot(t40$low_qual, t40$high_qual, main=paste("n=",length(t40$CDR3)))

a=100 * length(which(t10$related_to_top100 >= 1)) / length(t10$CDR3)
b=100 * length(which(t20$related_to_top100 >= 1)) / length(t20$CDR3)
c=100 * length(which(t30$related_to_top100 >= 1)) / length(t30$CDR3)
d=100 * length(which(t40$related_to_top100 >= 1)) / length(t40$CDR3)

barplot(c(a,b,c,d), names=c("10", "20", "30", "40"), xlab="Base quality threshold", ylab = "Has similarity to top 100 cdr3's (%)")
