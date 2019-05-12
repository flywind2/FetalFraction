library(data.table) # for fast fread
library(dplyr)

normalize <- function(df){
  df <- data.frame(df)
  m <- as.matrix(df[,-1])
  rownames(m) <- df[,1]    
  #return(m[,5:28])
  m <- sweep(m[,5:28],1,m[,1],`/`)
  return(m)
}

header <- colnames(fread("headery.txt", header=T,sep = ";"))
panel <- fread("bins_questions_50000_q35_all.csv", header=F)
colnames(panel) <- header 

males <- read.csv2("males.txt")
females <- read.csv2("females.txt")

panel$patient <- gsub('_-', '-', panel$patient)
panel$patient <- gsub('_', '-', panel$patient)
panel$patient <- gsub('-[0-9][0-9][0-9][0-9]+', '', panel$patient)


# Keep male fetuses
#m <- filter(panel, panel$patient %in% males$patient)
#f <- filter(panel, panel$patient %in% females$patient)

males_panel <- colSums(normalize(filter(panel, panel$patient %in% males$patient)))
females_panel <- colSums(normalize(filter(panel, panel$patient %in% females$patient)))
newdf <- rbind(as.data.frame(t(males_panel), row.names = "males"), as.data.frame(t(females_panel),row.names = "females"))

#write.csv2(newdf,"excel")

barplot(as.matrix(newdf),cex.names=0.5,main="Read count distribution on Y",las=2,beside=T,legend = rownames(newdf),col=c("black","white"))
barplot(as.matrix(data.frame(newdf$X7275959.7562959,newdf$X6857959.7147959)),cex.names=1,main="Read count distribution on 7275959 and 6857959",las=1,beside=T,legend = rownames(newdf),col=c("black","white"))
