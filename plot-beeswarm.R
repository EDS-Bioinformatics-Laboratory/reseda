library('beeswarm')

myfile="/home/barbera/git/lineage-tree/beeswarm10.png.csv"
d<-read.csv(myfile, sep=" ", header = FALSE)
colnames(d)<-c("seq","freq","frac","color","accs")
d$sample="D10"
attach(d)

# Enlarge margins: mar=c(bottom, left, top, right)
# Allow to draw outside the plot region: xpd=TRUE
# vertical axis labels: las=2
par(mar=c(8,4,4,2), xpd=TRUE)

# Make the beeswarm plot
beeswarm(frac ~ sample, data=d,
         log = TRUE, pch = 16,
         pwcol=color,
         main = 'Subclones', method="swarm", corral="wrap", corralWidth=0.7, # was 0.9
         las=2)
# Add legend
coord<-locator(1)
#legend(coord$x,coord$y, legend = c("any","GEM1","MAIT","GEM2"),
#       title = "Cell types", pch = 16, col = c("#00000050",2,3,4),cex=0.75)
legend(coord$x,coord$y, legend = c("GEM1","MAIT","GEM2"),
       title = "Cell types", pch = 16, col = c(2,3,4),cex=0.75)

par(def.par)#- reset to default
