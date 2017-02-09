# New contamination script for Miseq

# read.me for PLK
### now it is a loop for top1,10,25, and 100 clones
### requires clones file
### IT IS CRITICAL THAT CLONES$ID2 HAS PT [space] SAMPLE

# input
clones.file = "/home/barbera/Data/tbcell/lineage-trees/20161201-SP-RAII-rituximab/2016-12-01-RA2_BCRh_PB_Tree-clones.csv"

# Read clones file and add the ID2 column: PT [space] SAMPLE
clones = read.csv(clones.file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
project = "RAII"

clones$ID2 = paste(clones$pt,clones$sample)
clones$ID2=as.factor(clones$ID2)

clones$VJCDR3 = paste(paste(clones$V.gene, clones$J.gene, sep = "-"), clones$cdr3pep, sep = "-")

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
