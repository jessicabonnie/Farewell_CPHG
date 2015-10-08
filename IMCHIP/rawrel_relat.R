print("#########################################################################")
print("#")
print("# GWAS-Pipedream - Standard GWAS QC Pipeline Package")
print("# (c) 2014-2015 JBonnie, WMChen")
print("# Script: rawrel_relat.R -- colors by relationship")
print("# Usage: rawrel_relat.R <kinship table path> <nickname for prefixing file> <Overall Study Title> <2nd Degree Kinship Threshold>")
print("#")
print("#########################################################################")

args = commandArgs(TRUE)
print(args)
table_loc <- as.character(args[1])
nickname <- as.character(args[2])
overall_title <- gsub("_"," ",as.character(args[3]))


dupthresh <- 0.4
PO_IBS <- .005
FSthresh <- 0.18
SDthresh <-0.1


data <- read.table(file = table_loc, header=T)
N <- sum(data$Kinship > SDthresh)


postscript(paste0(nickname,"rawrel_relat.ps"), paper="letter", horizontal=T)
#jpeg(paste0(nickname,"rawrel_relat.jpg"),width=1000, height=1000)
plot(data$IBS0, data$Kinship, type="p",
     col = "lightgray", cex.lab=1.3,cex=1.5,
     #main = paste0("Unreported Relationships in ", overall_title,"\n(",N," Pairs with Kinship > 0.0884)"),
     main = paste0("Relationships in ", overall_title,"\n(",N," Pairs with Kinship > ",as.character(SDthresh),")"),
     xlab="Proportion of Zero IBS", ylab = "Estimated Kinship Coefficient", xlim=c(0,max(data$IBS0)),ylim=c(0,0.5))


points(data$IBS0[data$Kinship > dupthresh],data$Kinship[data$Kinship > dupthresh],
       col="red",cex=1.5)
points(data$IBS0[data$IBS0 <= PO_IBS & data$Kinship < dupthresh],data$Kinship[data$IBS0 <= PO_IBS & data$Kinship < dupthresh],
       col="darkgoldenrod1",cex=1.5)
points(data$IBS0[data$IBS0 > PO_IBS & data$Kinship >= FSthresh],data$Kinship[data$IBS0 > PO_IBS & data$Kinship >= FSthresh],
       col="chartreuse3",cex=1.5)
points(data$IBS0[data$Kinship < FSthresh & data$Kinship >= SDthresh],data$Kinship[data$Kinship < FSthresh & data$Kinship >= SDthresh],
       col="blue",cex=1.5)

u <- par("usr")
legend(u[2], u[4], xjust=1, yjust=1, 
       legend=c("Duplicate / MZ Twin","Parent-Offspring","Full Sibling","2nd Degree"),
       col=c("red","darkgoldenrod1","chartreuse3","blue"),
       #text.col = c("magenta","lightgray"),
       pch = "O", cex = 1.5)
#pch = 19, cex = 1.2)

dev.off()


