
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



data <- read.table(file=datafile, header=F, na.string="X")
dn <- data$V6==2

postscript(file="EURprojpca.ps", paper="letter", horizontal=T)
plot(data[,7], data[,8], type="p", xlab="PC1", ylab="PC2", main = paste("CEU/TSI + ",overall_title), col = "black")
points(data[dn,7], data[dn,8], col = "red")
abline(v=0, col="green")

u <- par("usr")
legend(u[2], u[4], xjust=1, yjust=1, 
legend=c("CEU+TSI", overall_title),  col=c("black","red"),text.col = c("black", "red"), pch = 19, cex = 1.2)

dev.off()
