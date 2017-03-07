# RESEDA - REPertoire SEquencing Data Analysis
# Copyright (C) 2016 Paul L Klarenbeek
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Contamination script for TCR/BCR repertoire data

### Compare top1, 10, 25, and 100 clones between all samples
### Input: clones.csv file
### Requirements: clones$ID2 should contain:  PT [space] SAMPLE (patient sample-of-patient)
###               clones$VJCDR3 contains the clone name, e.g.: V-J-CDR3
### Output: <date>.<projectname>.overlap.plots.pdf

# Input file
clones.file1 = "/mnt/immunogenomics/RUNS/run13-20170224-miseq/results-tbcell/run-clones_subs-VDJmouse-IGK_MOUSE.csv"
clones.file2 = "/mnt/immunogenomics/RUNS/run13-20170224-miseq/results-tbcell/run-clones_subs-VDJmouse-TRB_MOUSE.csv"
clones.file3 = "/mnt/immunogenomics/RUNS/run13-20170224-miseq/results-tbcell/run-clones_subs-C1IgG4-IGH_HUMAN.csv"
clones.file4 = "/mnt/immunogenomics/RUNS/run13-20170224-miseq/results-tbcell/run-clones_subs-Paired-IGH_HUMAN.csv"
pt.file = "Run013-pt-table.csv"

# Read clones file
clones1 = read.csv(clones.file1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
clones2 = read.csv(clones.file2, header=TRUE, sep="\t", stringsAsFactors=FALSE)
clones3 = read.csv(clones.file3, header=TRUE, sep="\t", stringsAsFactors=FALSE)
clones4 = read.csv(clones.file4, header=TRUE, sep="\t", stringsAsFactors=FALSE)
pt.table = read.csv(pt.file, header=TRUE, sep=",", stringsAsFactors=FALSE)

clones = rbind(clones3, clones4)
#clones = clones1
clones$SampleName=gsub("_S[0-9]+","",clones$Sample)
clones = merge(clones,pt.table, by.x="SampleName", by.y="ID", all.x = TRUE)
project = "run13-IGH_HUMAN"

# Add the ID2 column: PT [space] SAMPLE
clones$ID2 = paste(clones$pt,clones$SampleName)
#clones$ID2 = clones$Sample
clones$ID2=as.factor(clones$ID2)


# Add VJCDR3 column if it is not already present: V-J-CDR3
#clones$VJCDR3 = paste(clones$V.gene, clones$J.gene, clones$cdr3pep, sep = "-")
clones$VJCDR3 = paste(clones$V_sub, clones$J_sub, clones$cdr3pep, sep = "-")


########## MAIN ############

########## make name
name=paste(Sys.Date(),project,'overlap.plots','pdf',sep='.')

########## start pdf
pdf(name)


# loop
numbers=c(1,10,25,100)
for (ii in 1:4){
    start=Sys.time()
    n=numbers[ii]
    title=paste('overlap of', n, 'most dominant clones')


    # create clones2 (only top n of each sample)
    output=vector('list')                                  # output vector
    ID=unique(clones$ID2)                                  # unique IDs
    for(i in 1:length(ID)){                                # loop
        temp1=clones[clones$ID2==ID[i],]                   # subset data per ID
        temp2=temp1[1:n,]                                  # select n clones
        output[[i]]=temp2                                  # collect in vector, end loop
    }
    clones2=do.call('rbind',output)                        # make matrix from output (clones2)
    loop=as.character(sort(unique(clones2$ID2)))           # loop of ID2
    clones2$VJCDR3=paste(clones2$VJCDR3,'x')               # add x to prevent grep problems

    # overlap with all samples
    output=vector('list')                                       # output vector
    for (i in 1:length(loop)){                                  # loop
        temp=clones2$VJCDR3[as.character(clones2$ID2)==loop[i]] # select clones from 1 sample,
        temp2=clones2$ID2[clones2$VJCDR3 %in% temp]
        temp3=table(temp2)                                      # table
        temp3[temp3>n] <- n                                     # substitute samples >n with n
        output[[i]]=temp3
    }
    output2=do.call('rbind',output)
    rownames(output2)=loop

    # figure
    ### palette and borders
    if(n==1){
        my_palette <- colorRampPalette(c("white","blue"))(n = 2)
        col_breaks = c(0,0.99,1)
    }
    if(n>=10){
        my_palette <- colorRampPalette(c("white","blue"))(n = n)
        col_breaks = seq(0,n,1)
    }

    ### make plot
    image(output2,
          col=my_palette,
          breaks=col_breaks,axes=F,
          main=title)
    ### boxes
    temp1=gsub(' .*','',rownames(output2))                   # adjust rownames
    temp2=unique(temp1)                                      # extract unique pt's
    temp3=match(temp2,temp1)                                 # find positions of unique pts
    temp4=seq(0,1,by=1/(nrow(output2)-1))                    # find sample positions on X-axis (0:1)
    temp5=temp4[temp3]                                       # find middle of box start positions
    temp6=temp5-((1/nrow(output2))/2)                        # find left border of box start positions
    temp7=c(temp6[2:length(temp6)],1+((1/nrow(output2))/2))  # find right border of box start positions
    ### draw boxes
    for (i in 1:length(temp6)){
        x0=temp6[i]; x1=temp7[i];
        y0=x0; y1=x1
        segments(x0=x0,x1=x1,y0=y0,y1=y0)
        segments(x0=x0,x1=x1,y0=y1,y1=y1)
        segments(x0=x0,x1=x0,y0=y0,y1=y1)
        segments(x0=x1,x1=x1,y0=y0,y1=y1)
    }
    ### axes
    axis(2,at=temp5,labels=temp2,las=2, cex.axis=0.5)
    axis(1,at=temp5,labels=temp2,las=2, cex.axis=0.5)

    # print time taken
    end=Sys.time()
    print(paste('figure top', n, 'clones done in', round((end-start), digits=0), 'seconds'))

    # end loop
}

# close pdf
dev.off()

cat(name, "written to disk\n")
