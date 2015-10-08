
args = commandArgs(TRUE)
print(args)
datafile <- as.character(args[1])
overall_title <- gsub("_"," ",as.character(args[2]))
chip <- as.character(args[3])



if (chip == "I"){ 
  chipname <- "Immunochip"
}else if (chip == "E"){ 
  chipname <- "HumanCoreExomeChip"
}else if (chip %in% c("12v1.0","12v1.1","24v1.0","24v1.1","24v1","12v1")){ 
  chipname <- paste("HumanCoreExomeChip",chip)}

if (is.null(chip)){ chipname <- ""}


        print(datafile)


data <- read.table(file = datafile, header=T)

postscript("samplecallrate.ps", paper="letter", horizontal=T)
par(mfrow=c(2,1))
x <- 1-data$Missing
hist(x, xlab = "Call Rate", main = paste("Histogram of Sample Call Rate in", overall_title,chipname, "Data"))
hist(x[x>=0.992 & x < 0.9995], xlim=c(0.992,0.9995), xlab = "Call Rate", main = "Call Rate Between 99.2% and 99.95%")

par(mfrow=c(1,1))
x <- data$Heterozygosity
hist(x, xlab="Heterozygosity", main=paste("Histogram of Heterozygosity in", overall_title,chipname, "Data"))

x <- 1-data$Missing
y <- data$Heterozygosity
plot(x, y, xlab="Call Rate", ylab="Heterozygosity", main=paste("Sample Level QC in", overall_title,chipname, "Data"), col="green")
abline(v=0.95, col="red")
dev.off()

