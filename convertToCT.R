install.packages("protr");
library("protr")
library("readr")

PPIData <- read_delim("Downloads/11050", ";", escape_double = FALSE, trim_ws = TRUE)

mentha = readFASTA(system.file(args[1], package = "protr"))
names <- substr(names(mentha), 4,9)

toStringFromVector <- function(v,w,x) {
  str <- x
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
  if((PPIData$`Protein A`[i] %in% names) && (PPIData$`Protein B`[i] %in% names)){
    ProteinA <- extractCTriad(mentha[[getIndex(PPIData$`Protein A`[i])]])
    ProteinB <- extractCTriad(mentha[[getIndex(PPIData$`Protein B`[i])]])
    vectorOfData[i] <- toStringFromVector(ProteinA, ProteinB, "1")
  }
  else{
  }
}

vectorOfNegativeData <- vector(mode="character", length=nrow(PPIData))

for(i in 1:nrow(PPIData)){
  index <- sample(1:(nrow(PPIData)),1)
  print(index)
  index2 <- sample(1:(nrow(PPIData)), 1)
  print(index2)
  
  if((PPIData$`Protein A`[index] %in% names) && (PPIData$`Protein B`[index2] %in% names)){
    print("in")
    ProteinA <- extractCTriad(mentha[[getIndex(PPIData$`Protein A`[index])]])
    ProteinB <- extractCTriad(mentha[[getIndex(PPIData$`Protein B`[index2])]])
    vectorOfNegativeData[i] <- toStringFromVector(ProteinA, ProteinB, "-1")
  }
  else{
  }
}


fileConn<-file("~/Documents/GitHub/zikaHumanPPI/data/positivePPI.txt")
writeLines(vectorOfData, fileConn)
close(fileConn)

fileConn<-file("~/Documents/GitHub/zikaHumanPPI/data/negativePPI.txt")
writeLines(vectorOfNegativeData, fileConn)
close(fileConn)

