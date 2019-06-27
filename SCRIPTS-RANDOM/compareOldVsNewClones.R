library(plyr)

setwd("git/tbcell-miseq-pipeline/")

# Note: the clone definition here is: V+J+CDR3pep

# Function for the mode. From: https://www.r-bloggers.com/computing-the-mode-in-r/
Mode = function(x){
  ta = table(x)
  tam = max(ta)
  if (all(ta == tam))
    mod = NA
  else
    if(is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
  else
    mod = names(ta)[ta == tam]
  return(mod)
}

# Read allinfo file with mutation information, filter it and create a clones dataframe
df.allinfo = read.csv("final/SP-CB21-Bu_S176-AGCTAGCT-mutations-per-accession.csv", sep=",", header = T, stringsAsFactors = F)
df.allinfo.filtered=df.allinfo[which(df.allinfo$V_sub!='None' & df.allinfo$J_sub!='None' & df.allinfo$cdr3_qual_min>=30 & (df.allinfo$V_flag == 0 | df.allinfo$V_flag == 16) & (df.allinfo$J_flag == 0 | df.allinfo$J_flag == 16)),]
clones.old=ddply(df.allinfo.filtered, .(V_sub,J_sub,cdr3pep), summarise, freq=length(unique(acc)), umis=length(unique(beforeMID)), total.mut.x=sum(mut.count.x), median.mut.x=median(mut.count.x), avg.mut.x=mean(mut.count.x), mode.mut.x=min(Mode(mut.count.x)), avg.mutfrac.x=mean(mut.frac.x), align.avg=mean(align.length.x), align.min=min(align.length.x), align.max=max(align.length.x), align.std=sd(align.length.x))

# Read the new clones file
clones.new=read.csv("final/SP-CB21-Bu_S176_L001.assembled-AGCTAGCT-IGH_HUMAN-clones-mut-sites-reassigned.csv", sep='\t', header=T, stringsAsFactors = F)

# Merge into one dataframe
df=merge(clones.old, clones.new, by.x = c('cdr3pep','V_sub','J_sub'), by.y = c('cdr3pep','V_sub','J_sub'))

# Compare Frequency, UMIs, total mutation count and mutation fraction
pdf("SampleGraph.pdf",width=7,height=7)
plot(df$acc.nunique, df$freq, log='xy', main="NEW vs OLD frequency")
plot(df$beforeMID.nunique, df$umis, log='xy', main="NEW vs OLD UMIs")
plot(df$mut.count_x.sum, df$total.mut.x, log='xy', main="NEW vs OLD sum of mutation count")
plot(df$mut.frac_x.mean, df$avg.mutfrac.x, main="NEW vs OLD average mutation fraction\naverage(nr of mutations / alignment length)")

# Compare mode of mutations with the total/freq, the median and the average
plot(df$median.mut.x, df$avg.mut.x, main="OLD median vs average")
plot(df$mode.mut.x, df$avg.mut.x, main="OLD mode vs average")
plot(df$mode.mut.x, df$median.mut.x, main="OLD mode vs median")
plot(df$mut.count_x..lambda., df$avg.mut.x, main="NEW mode vs OLD average")
plot(df$mut.count_x..lambda., df$median.mut.x, main="NEW mode vs OLD median")
plot(df$mut.count_x..lambda., df$mode.mut.x, main="NEW vs OLD mode")
dev.off()

# Check why the average fraction (mut/alignmentlength) is so different
boxplot(clones.old$align.avg, clones.old$align.min, clones.old$align.max, names=c('align.avg', 'align.min', 'align.max'))
test=clones.old$align.std
test[is.na(test)]=0
boxplot(test)

# Check if the alignment length differs a lot with same VJCDR3+UMI. Answer: YES
vjcdr3umi=ddply(df.allinfo.filtered, .(V_sub,J_sub,cdr3pep,beforeMID), summarise,freq=length(unique(acc)), align.count=length(unique(align.length.x)), mut.count=sum(mut.count.x), align.avg=mean(align.length.x), align.min=min(align.length.x), align.max=max(align.length.x), align.std=sd(align.length.x))
test2=vjcdr3umi$align.std
test2[is.na(test2)]=0
boxplot(test2)
boxplot(vjcdr3umi$align.std, main="VJCDR3-UMI standard deviation of alignment length")

# Check if number of mutations goes up if the alignment length is longer (that should be the case). Answer: YES
plot(df.allinfo.filtered$align.length.x, df.allinfo.filtered$mut.count.x)

# Check if fraction of mutations is about the same if the alignment length gets longer. Answer: plot looks biased with extremes for shorter reads
plot(df.allinfo.filtered$align.length.x, df.allinfo.filtered$mut.frac.x)
