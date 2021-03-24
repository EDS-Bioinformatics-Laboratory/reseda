# Example from commandline:
# Rscript PlotVJ.R myfile=\"BASTA-2202-B_S57_L001.assembled-ACGTACGT-IGH_HUMAN-clones-mut-sites-reassigned.csv\"

# Batch processing:
# myfiles=`ls *.csv`; for myfile in $myfiles; do Rscript PlotVJ.R myfile=\"$myfile\"; done
# mogrify -density 96 -flatten -format png *.pdf

args = (commandArgs(TRUE))

# Input parameters
if(length(args)==0){
    cat("ERROR: no arguments supplied\n")
    cat("  myfile='SAMPLE_Sn_L001-etc-clones-mut-sites-reassigned.csv'\n")
    q()
} else {
    for(i in 1:length(args)){
        cat(i, args[[i]], "\n")
        eval(parse(text=args[[i]]))
    }
    cat("arguments given to script\n")
    # q()
}

library(circlize)
library(reshape)

# Read the file
mytitle = sub("_L001.*", "", myfile)
df = read.csv(myfile, sep="\t", header=T, stringsAsFactors = F)

# Select columns
Vgene=df$V_sub
Jgene=df$J_sub
UMIs.frac=df$UMIs.frac  # NOTE: if there are no UMIs you should use the df$freq column instead

# Select only the first V and J assignment
Vgene = sub("\\+.*", "", Vgene)
Jgene = sub(",.*", "", Jgene)

# Reformat data
mat <- data.frame(Vgene,Jgene,UMIs.frac)
#mat <- with(mat, table(Vgene, Jgene))
mat<-cast(mat, Vgene ~ Jgene, fun.aggregate=sum)
rownames(mat) <- mat[,1]
mat[,1] <- NULL
mat<-as.matrix(mat)

# Make the circular plot
pdf(paste(mytitle, "circos.pdf", sep="-"),width=7,height=7)
grid.col <- setNames(rainbow(length(unlist(dimnames(mat)))), union(rownames(mat), colnames(mat)))
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)

# now, the image with rotated labels
chordDiagram(mat, annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)
title(main=mytitle)
dev.off()

#circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
