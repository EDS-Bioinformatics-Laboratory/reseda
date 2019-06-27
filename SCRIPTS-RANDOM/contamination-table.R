#### script for contamination table

# read.me
# used for preRA.TCR, might need adjustment in the future
# I made a table of the top 10 overlapping clones, which was enough to spot contamination problems
# the script finds all clones that are in > 1 patient
# it produces output4 which is a list with the overlapping clones

n=10
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

######## top10 select
# select clones that are in >1 patient
output2=vector(mode='numeric')
for(i in 1:length(clones2$VJCDR3)){
    temp=clones2$VJCDR3[i]
    output[i]=length(unique(clones2$pt[clones2$VJCDR3==temp]))
    print(i)
}
temp2=output>1
temp3=grep('TRUE',temp2)
temp4=unique(clones2$VJCDR3[temp3])

# find corresponding rows
output=vector('list')
for (i in 1:length(temp4)){
    output[[i]]=clones2[clones2$VJCDR3==temp4[i],]
}
output2=do.call('rbind',output)

# organize sample per clone
output3=vector('list')
loop=unique(output2$VJCDR3)
for (i in 1:length(loop)){
    temp=output2[output2$VJCDR3==loop[i],]
    temp2=paste(unique(paste(temp$pt,temp$sample,temp$ID,sep='.')),collapse='__')
    output3[[i]]=c(temp$VJCDR3[1],temp2)
}
output4=do.call('rbind',output3)



n=25
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

######## top25 select
# select clones that are in >1 patient
output2=vector(mode='numeric')
for(i in 1:length(clones2$VJCDR3)){
    temp=clones2$VJCDR3[i]
    output[i]=length(unique(clones2$pt[clones2$VJCDR3==temp]))
    print(i)
}
temp2=output>1
temp3=grep('TRUE',temp2)
temp4=unique(clones2$VJCDR3[temp3])

# find corresponding rows
output=vector('list')
for (i in 1:length(temp4)){
    output[[i]]=clones2[clones2$VJCDR3==temp4[i],]
}
output2=do.call('rbind',output)

# organize sample per clone
output3=vector('list')
loop=unique(output2$VJCDR3)
for (i in 1:length(loop)){
    temp=output2[output2$VJCDR3==loop[i],]
    temp2=paste(unique(paste(temp$pt,temp$sample,temp$ID,sep='.')),collapse='__')
    output3[[i]]=c(temp$VJCDR3[1],temp2)
}
output4=do.call('rbind',output3)
