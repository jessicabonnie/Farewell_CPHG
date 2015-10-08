
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

postscript("snpqc.ps", paper="letter", horizontal=T)
par(mfrow=c(3,1))
x <- data$CallRate
range.plot <- x>0.95 & x < 0.995 & data$Chr!="Y"
#hist(x, xlab = "Call Rate", main = "Histogram of SNP Call Rate in JDRF_DN ExomeChip Data")
hist(x[range.plot], xlim=c(0.95,0.995), xlab = "Call Rate", main = "Call Rate Between 95% and 99.5%")     
variant.notcommon <- data$Freq_A < 0.05 | data$Freq_A > 0.95
hist(x[range.plot & variant.notcommon], xlim=c(0.95,0.995), xlab = "Call Rate", main = "Call Rate Between 95% and 99.5% for SNPs with MAF < 0.05")
variant.rare <- data$Freq_A < 0.01 | data$Freq_A > 0.99
hist(x[range.plot & variant.rare], xlim=c(0.95,0.995), xlab = "Call Rate", main = "Call Rate Between 95% and 99.5% for SNPs with MAF < 0.01")

par(mfrow=c(2,1))
x <- 0.5 - abs(data$Freq_A-0.5)
hist(x, xlab = "MAF", main = paste("Histogram of MAF in", overall_title,chipname, "Data"))
hist(x[x<0.05], xlab = "MAF", main = "Histogram of MAF in Rare Variants")

par(mfrow=c(1,1))
N.MZ <- max(data$N_MZ)
x <- data$Err_InHetMZ
y <- data$Err_InMZ
variant.error <- data$N_errMZ > 0
variant.rare <- data$Freq_A < 0.01 | data$Freq_A > 0.99
variant.lessfreq <- (data$Freq_A < 0.05 & data$Freq_A >= 0.01) | (data$Freq_A <= 0.99 & data$Freq_A > 0.95)
plot(x[variant.error], y[variant.error], xlab = "Adjusted Error Rate", ylab="Error Rate", main=paste("Problematic SNPs In", N.MZ, "Pairs of Duplicates"), col="green")
points(x[variant.error & variant.lessfreq], y[variant.error & variant.lessfreq], col="blue")
points(x[variant.error & variant.rare], y[variant.error & variant.rare],  col="red")
abline(v=0.1, col="purple")
abline(h=0.01, col="purple")
u <- par("usr")
legend(u[2], u[4], xjust=1, yjust=1,
legend=c("Common Variants", "Less Frequent Variants", "Rare Variants"), col=c("green", "blue", "red"),  
text.col = c("green", "blue", "red"), pch = 19, cex = 1.2) 

x <- data$CallRate
y <- data$Err_InHetMZ
variant.error <- (data$N_errMZ > 0 & data$Chr!="Y")
plot(x[variant.error], y[variant.error], xlab = "Call Rate", ylab = "Adjusted Error Rate", main=paste("Problematic SNPs In", N.MZ, "Pairs of Duplicates"), col="green")
points(x[variant.error & variant.lessfreq], y[variant.error & variant.lessfreq], col="blue")
points(x[variant.error & variant.rare], y[variant.error & variant.rare],  col="red")
abline(v=0.95, col="green")
abline(v=0.99, col="red")
abline(h=0.1, col="purple")
u <- par("usr")
legend(u[1], u[4], xjust=0, yjust=1,
legend=c("Common Variants", "Less Frequent Variants", "Rare Variants"), col=c("green", "blue", "red"),
text.col = c("green", "blue", "red"), pch = 19, cex = 1.2)

N.PO <- max(data$N_PO)
x <- data$Err_InHomPO
y <- data$Err_InPO
variant.error <- data$N_errPO > 0
variant.rare <- data$Freq_A < 0.01 | data$Freq_A > 0.99
variant.lessfreq <- (data$Freq_A < 0.05 & data$Freq_A >= 0.01) | (data$Freq_A <= 0.99 & data$Freq_A > 0.95)
plot(x[variant.error], y[variant.error], xlab = "Adjusted Error Rate", ylab="Error Rate", main=paste("Problematic SNPs In", N.PO, "Parent-Offspring Pairs"), col="green")
points(x[variant.error & variant.lessfreq], y[variant.error & variant.lessfreq], col="blue")
points(x[variant.error & variant.rare], y[variant.error & variant.rare],  col="red")
abline(v=0.1, col="purple")
abline(h=0.01, col="purple")
u <- par("usr")
legend(u[2], u[3], xjust=1, yjust=0,
legend=c("Common Variants", "Less Frequent Variants", "Rare Variants"), col=c("green", "blue", "red"),
text.col = c("green", "blue", "red"), pch = 19, cex = 1.2)

N.Trio <- max(data$N_trio)
x <- data$Err_InHetTrio
y <- data$Err_InTrio
variant.error <- data$N_errTrio > 0
variant.rare <- data$Freq_A < 0.01 | data$Freq_A > 0.99
variant.lessfreq <- (data$Freq_A < 0.05 & data$Freq_A >= 0.01) | (data$Freq_A <= 0.99 & data$Freq_A > 0.95)
plot(x[variant.error], y[variant.error], xlab = "Adjusted Error Rate", ylab="Error Rate", main=paste("Problematic SNPs In", N.Trio, "Parents-Offspring Trios"), col="green")
points(x[variant.error & variant.lessfreq], y[variant.error & variant.lessfreq], col="blue")
points(x[variant.error & variant.rare], y[variant.error & variant.rare],  col="red")
abline(v=0.1, col="purple")
abline(h=0.01, col="purple")
u <- par("usr")
legend(u[2], u[4], xjust=1, yjust=1,
legend=c("Common Variants", "Less Frequent Variants", "Rare Variants"), col=c("green", "blue", "red"),
text.col = c("green", "blue", "red"), pch = 19, cex = 1.2)

dev.off()

