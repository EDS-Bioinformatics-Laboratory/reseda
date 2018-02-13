library(httr)

# N.B. Create a variable 'passwd' with the password manually

url="https://researchdrive.surfsara.nl/remote.php/webdav/amc-immunogenomics/RUNS/run20-20171127-miseq/results-tbcell20171129/report-all.csv"
fh=GET(url, authenticate("bschaik", passwd))

if (http_status(fh)$reason == "OK") {
  df <- read.csv(textConnection(content(fh, 'text')), header=T, sep=" ")  
}
