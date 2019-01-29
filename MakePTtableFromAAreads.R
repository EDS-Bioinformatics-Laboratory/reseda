library("plyr")
aa.reads=read.csv("/mnt/immunogenomics/RUNS/runNN-20171122-RA2-ST-BCRh/2011-11-15-RA2-project.66.84.88.98-BCRh-reads.csv", header=T, sep=",", stringsAsFactors = F)
pt.table=ddply(aa.reads, .(REGIONMID.RUN,ID,PROJECT,pt), summarise, freq=length(unique(accession)))
pt.table$group=gsub(",", "-", pt.table$REGIONMID.RUN)
write.csv(pt.table, file="pt.table.20171122.RA2-ST-BCRh.csv")

aa.reads=read.csv("/mnt/immunogenomics/RUNS/runNN-20171122-RA2-ST-TCRb/2015-02-26_AA.reads-TCRb-RA2.csv", header=T, sep=",", stringsAsFactors = F)
pt.table=ddply(aa.reads, .(REGIONMID.RUN,MID,ID,project,pt,sample,group), summarise, freq=length(unique(accession)))
write.csv(pt.table, file="pt.table.20171122.RA2-ST-TCRb.csv")
