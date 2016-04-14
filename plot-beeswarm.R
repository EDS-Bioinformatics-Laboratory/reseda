library('beeswarm')

# VOOR ANTOINE: B cell lineages
#myfile="/home/barbera/git/lineage-tree/voor-antoine-b-cell-lineage/LN25-top3-lineages/V3.74-J4/beeswarm2.png.csv"
#d<-read.table(myfile, sep=" ", header = FALSE, stringsAsFactors=FALSE)
#colnames(d)<-c("seq","freq","frac","color","accs")
#d$sample="D2"
#attach(d)

# VOOR PIPELINE PAPER: clone figuren voor en na V re-assignment
d1<-read.table("/home/barbera/TMP/S074-002_S73_L001.assembled-ACGTACGT-TRB_HUMAN-clones-subs.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
d1$sample<-"S73-before"
d2<-read.table("/home/barbera/TMP/S074-002_S73_L001.assembled-ACGTACGT-TRB_HUMAN-all_info.csv.rr.clones_subs.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
d2$sample<-"S73-after"
d3<-read.table("/home/barbera/TMP/S074-004_S74_L001.assembled-ACTGACTG-TRB_HUMAN-clones-subs.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
d3$sample<-"S74-before"
d4<-read.table("/home/barbera/TMP/S074-004_S74_L001.assembled-ACTGACTG-TRB_HUMAN-all_info.csv.rr.clones_subs.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
d4$sample<-"S74-after"

# Concatenate everything
d<-rbind(d3,d4)  # bind everything
d$color="black"
d$color[which(d$V_sub=="TRBV12-3")]="red"
d$color[which(d$V_sub=="TRBV12-4")]="red"
d$color[which(d$V_sub=="TRBV12-3+TRBV12-4")]="green"

# Remove 90% of entries with read frequency == 1
inx_1<-which(d$freq == 1)
s_inx<-sample(inx, 0.1*length(inx))
inx_higher<-which(d$freq > 1)
keep<-sort(c(inx_higher,s_inx))
d<-d[keep,]

# Enlarge margins: mar=c(bottom, left, top, right)
# Allow to draw outside the plot region: xpd=TRUE
# vertical axis labels: las=2
#par(mar=c(8,4,4,2), xpd=TRUE)

# Make the beeswarm plot
attach(d)
beeswarm(freq ~ sample, data=d,
         log = TRUE, pch = 16,
         pwcol=color,
         main = 'Clones before and after V gene re-assignment', method="swarm", corral="wrap", corralWidth=0.7, # was 0.9
         las=2)

# Add legend
#coord<-locator(1)
#legend(coord$x,coord$y, legend = c("any","GEM1","MAIT","GEM2"),
#       title = "Cell types", pch = 16, col = c("#00000050",2,3,4),cex=0.75)
#legend(coord$x,coord$y, legend = c("GEM1","MAIT","GEM2"),
#       title = "Cell types", pch = 16, col = c(2,3,4),cex=0.75)

#par(def.par)#- reset to default
