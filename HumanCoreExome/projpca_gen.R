print("#########################################################################")
print("#")
print("# GWAS-Pipedream - Standard GWAS QC Pipeline Package")
print("# (c) 2014-2015 JBonnie, WMChen")
print("# Script: projpca_gen.R")
print("# Usage: projpca_gen.R <pca-projection file> <Overall Study Title> <chip> <boundary on pc2>")
print("#")
print("#########################################################################")

args = commandArgs(TRUE)
print(args)
datafile <- as.character(args[1])
overall_title <- gsub("_"," ",as.character(args[2]))
chip <- as.character(args[3])
pc2boundary <- 0.45
if (length(args >3)){
  pc2boundary <- as.numeric(args[4])
}


if (chip == "I"){ 
  chipname <- "Immunochip"
}else if (chip == "E"){ 
  chipname <- "HumanCoreExomeChip"
}else if (chip %in% c("12v1.0","12v1.1","24v1.0","24v1.1","24v1","12v1")){ 
  chipname <- paste("HumanCoreExomeChip",chip)}
if (is.null(chip)){ chipname <- ""}


#args=commandArgs()[7]

        # Extracting key user specified parameter from 'arg's
    #    datafile=substr(args, 2, nchar(args))

  #      print(datafile)


data <- read.table(file=datafile, header=F, na.string="X")
dn <- data$V6==2

postscript(file="projpca.ps", paper="letter", horizontal=T)
plot(data[,7], data[,8], type="p", xlab="PC1", ylab="PC2", main = paste("HapMapIII +", overall_title,"-",chipname), col = "black")
points(data[dn,7], data[dn,8], col = "red")
abline(h=pc2boundary, col="green")

u <- par("usr")
legend(u[2], u[4], xjust=1, yjust=1, 
legend=c("HapMap", overall_title),  col=c("black","red"),text.col = c("black", "red"), pch = 19, cex = 1.2)

dev.off()
