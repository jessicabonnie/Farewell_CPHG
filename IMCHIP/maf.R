print("#########################################################################")
print("#")
print("# GWAS-Pipedream - Standard GWAS QC Pipeline Package")
print("# (c) 2014-2015 JBonnie, WMChen")
print("#")
print("#########################################################################")


args = commandArgs(TRUE)
print(args)
table_loc=as.character(args[1])
graph_loc=as.character(args[2])
overall_title=as.character(args[3])
chip=as.character(args[4])

print(chip)
data <- read.table(file = table_loc, header=T)

postscript(graph_loc, paper="letter", horizontal=T)
par(mfrow=c(2,1))
x <- 0.5 - abs(data$Freq_A-0.5)

if (chip =="I"){
chip_title='ImmunoChip'
} else {
  chip_title='ExomeChip'
}

hist(x, xlab = "MAF", main = paste("Histogram of MAF in", overall_title, chip_title," Data"))

hist(x[x<0.02], xlab = "MAF", main = "Histogram of MAF in Rare Variants")
dev.off()

