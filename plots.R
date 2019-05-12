library(data.table) # For fast read
library(dplyr)

# Distribution of SeqFF fetal fraction estimates
all <- fread("ff_all.txt", header=T)
hist(all$FF,col="blue", breaks = 20,xlim=c(0,50), ylim=c(0,80))

# Comparison of the SeqFF and chromosome Y based fetal fraction estimates
comparison <- fread("y_vs_seqff.txt", header=T)
plot(comparison$FF, comparison$`Y based (CCHT)`, main="Y vs SeqFF", xlab="SeqFF", ylab="Y based", pch=19, col="blue")

# Processed data (GB) in each stage
data <- data.frame(Commodity = factor(c("Concatination step", "Alignment step", "Filtering (MAPQ 35)", "Chromosome Y based input"),
levels = c("Concatination step", "Alignment step", "Filtering (MAPQ 35)", "Chromosome Y based input")),
Production = c(501.2, 1700, 1100, 136))
barplot(data$Production, names = data$Commodity, main = "Data size GB", horiz=TRUE, xlim=c(0,2000), col="blue")

# Sum of read counts
raw_regions <- fread("sums.txt", header=T)
barplot(raw_regions$`Sum of raw read counts across all individuals`,  names = raw_regions$`Regions in chromosomal order` , main = "Sum of raw read counts", col="blue",las=2)

normalize <- function(df){
  df <- data.frame(df)
  m <- as.matrix(df[,-1])
  rownames(m) <- df[,1]    
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

# Keep male fetuses + normalize
males_panel <- colSums(normalize(filter(panel, panel$patient %in% males$patient)))
females_panel <- colSums(normalize(filter(panel, panel$patient %in% females$patient)))
newdf <- rbind(as.data.frame(t(males_panel), row.names = "males"), as.data.frame(t(females_panel),row.names = "females"))

# write.csv2(newdf,"excel")

# Read count distribution on XR7 and XR4
barplot(as.matrix(newdf),cex.names=0.5,main="Read count distribution on Y",las=2,beside=T,legend = rownames(newdf),col=c("black","white"))
barplot(as.matrix(data.frame(newdf$X7275959.7562959,newdf$X6857959.7147959)),cex.names=1,main="Read count distribution on XR7 and XR4",las=1,beside=T,legend = rownames(newdf),col=c("black","white"))

# Create histograms:
# Fetal fraction comparison in male and female pregnancies based on chromosome Y using x well covered unique regions
df <- data.frame(names = c("Males mean","Males median","Females mean","Females median"))
for (i in 1:24) { 
  from = paste("bins_50000_x",i,".csv",sep = "")
  
  header <- colnames(fread("header.50000.txt", header=T,sep = ";"))
  panel <- fread(from, header=F)
  panel <- panel[,c(1,5)]
  colnames(panel) <- header[c(1,5)]
  panel <- panel[panel$FF != 0,]
  
  males <- read.csv2("males.txt")
  females <- read.csv2("females.txt")
  
  panel$patient <- gsub('_-', '-', panel$patient)
  panel$patient <- gsub('_', '-', panel$patient)
  panel$patient <- gsub('-[0-9][0-9][0-9][0-9]+', '', panel$patient)
  
  males_panel <- filter(panel, panel$patient %in% males$patient)
  females_panel <- filter(panel, panel$patient %in% females$patient)
  
  df2 <- data.frame(from = c(mean(males_panel$FF),median(males_panel$FF),mean(females_panel$FF),median(females_panel$FF)))
  colnames(df2) <- from
  png(filename=paste(i,"png",sep = "."))
  hist(main = i,males_panel$FF,col="blue",breaks = 15,xlim=c(0,60), ylim = c(0,50))
  hist(females_panel$FF,col="pink",breaks = 15, add=T,xlim=c(0,60))
  dev.off()
  
  df <- cbind(df,df2)
}


