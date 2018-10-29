install.packages("protr");
library("protr")
library("readr")

PPIData <- read_delim("Downloads/11050", ";", escape_double = FALSE, trim_ws = TRUE)

mentha = readFASTA(system.file(args[1], package = "protr"))
names <- substr(names(mentha), 4,9)

toStringFromVector <- function(v,w) {
  str <- "1"
  for(i in 1:length(v)){
    str <- paste(str,paste0(i,":"))
    str <- paste0(str, v[i])
  }
  for(i in 1:length(w)){
    j <- length(v) + i
    str <- paste(str, paste0(j, ":"))
    str <- paste0(str, w[i])
  }
  return(str)
}

getIndex <- function(i) {
  return(which(names == i))
}

vectorOfData <- vector(mode="character", length=nrow(PPIData))

for(i in 1:nrow(PPIData)){
  if(getIndex(PPIData$`Protein A`[i]) != 0 && getIndex(PPIData$`Protein B`[i]) != 0){
    ProteinA <- extractCTriad(mentha[[getIndex(PPIData$`Protein A`[i])]])
    ProteinB <- extractCTriad(mentha[[getIndex(PPIData$`Protein B`[i])]])
    vectorOfData[i] <- toStringFromVector(ProteinA, ProteinB)
  }
  else{
    
  }
  
}



fileConn<-file("~\Documents\Github\ZikaHumanPPI\data\positivePPI.txt")
writeLines(vectorOfData, fileConn)
close(fileConn)


