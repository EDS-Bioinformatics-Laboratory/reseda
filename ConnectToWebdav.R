library(httr)

# N.B. Create the variables 'username' and 'passwd' with the password manually

readRemoteFile<-function(url, skip.nr){
  fh=GET(url, authenticate(username, passwd))
  
  if (http_status(fh)$reason == "OK") {
    df <- read.csv(textConnection(content(fh, 'text')), header=F, sep="\t", skip = skip.nr)  
  }
  return(df)
}

# url="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/run20-20171127-miseq/results-tbcell20171129/report-all.csv"
# mydata=readRemoteFile(url)

align.v.url="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/run20-20171127-miseq/results-tbcell20171129/For-Giulia/HD4-DNA-Jtail_S120_L001.assembled-ACGTACGT-IGHV_human-e-clean.sam"
align.j.url="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/run20-20171127-miseq/results-tbcell20171129/For-Giulia/HD4-DNA-Jtail_S120_L001.assembled-ACGTACGT-IGHJ_human-e-clean.sam"
align.v.again.url="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/run20-20171127-miseq/results-tbcell20171129/For-Giulia/HD4-DNA-Jtail_S120_L001.assembled-ACGTACGT-IGHV_human-unmapped-aligned-again.sam"
align.j.again.url="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/run20-20171127-miseq/results-tbcell20171129/For-Giulia/HD4-DNA-Jtail_S120_L001.assembled-ACGTACGT-IGHJ_human-unmapped-aligned-again.sam"

# Read data
d.v=readRemoteFile(align.v.url, 293)
d.j=readRemoteFile(align.j.url, 18)
d.v.again=readRemoteFile(align.v.again.url, 0)
d.j.again=readRemoteFile(align.j.again.url, 0)

# Get all original alignments with a V and J assigned
aligned=merge(d.v, d.j, by="V1")

# Combine V and J alignments (after second alignment round)
aligned.again=merge(d.v.again, d.j.again, by="V1", all.x=T)

# Get entries with a V, without a J
aligned.no.j=aligned.again[which(is.na(aligned.again$V2.y)),]

# Length of the sequences with and without J alignment
aligned$seq.length=nchar(as.vector(aligned$V10.x))
aligned.again$seq.length=nchar(as.vector(aligned.again$V10.x))
aligned.no.j$seq.length=nchar(as.vector(aligned.no.j$V10.x))
boxplot(aligned$seq.length, aligned.no.j$seq.length, names = c("Aligned VJ", "Aligned V again, no J"))

write.csv(aligned.no.j, file="HD4-DNA-Jtail-V-noJ.csv")
