print("#########################################################################")
print("#")
print("# GWAS-Pipedream - Standard GWAS QC Pipeline Package")
print("# (c) 2014-2015 JBonnie, WMChen")
print("#")
print("#########################################################################")


args = commandArgs(TRUE)
print(args)
table_loc=as.character(args[1])
hapmap_loc=as.character(args[2])
poplist_loc=as.character(args[3])
nickname=as.character(args[4])
overall_title=as.character(args[5])
covariablefile=as.character(args[6])
colorfile=as.character(args[7])
population=as.character(args[8])

covariables <- readLines(covariablefile)
colors <- readLines(colorfile)
covariable_num =length(covariables)
print(covariables)
#covariables <- read.table(file="../covariables.list")


draw_covariable <- function(covariable, covariables, data, colors){
  value <- as.numeric(which(covariables==covariable)[1])
  cv <-as.numeric(data$V9)==value
  points(data[cv,7], data[cv,8], col=colors[value])
  print(as.character(covariable))
  print(colors[value])
}

draw_covariable_sep <- function(covariable, covariables, data, data.hapmap, colors){
    plot(data.hapmap[,7], data.hapmap[,8], type="p", xlab="PC1", ylab="PC2", main = paste("HapMapIII +",overall_title,"-",covariable), col = "black")
    value <- as.numeric(which(covariables==covariable)[1])
    cv <-as.numeric(data$V9)==value
    #points(data[cv,7], data[cv,8], col=colors[value])
    points(data[cv,7], data[cv,8], col="purple")
  }


data.hapmap <- read.table(file=hapmap_loc, header=F, na.string="X")
data.pca <-read.table(file=table_loc, header=F, na.string="X")
poplist <- read.table(poplist_loc)
#data.hapmap <- read.table(file="../3_structure/hapmappc.txt", header=F, na.string="X")
#data.pca <-read.table(file="../3_structure/projpc.txt", header=F, na.string="X")
#nelist <- read.table("NElist.txt")
#data.pca[data.pca$V2 %in% nelist,]
postscript(file=paste0(nickname,"pca_",population,"_hwe.ps"), paper="letter", horizontal=T)
plot(data.hapmap[,7], data.hapmap[,8], type="p", xlab="PC1", ylab="PC2", main = paste0("HapMapIII + ",overall_title,", ",population), col = "black")

lapply(covariables, draw_covariable,covariables,data.pca[data.pca$V2 %in% poplist$V2,], colors)
u <- par("usr")


legend(u[1], u[4], 
       legend=c("HapMap",covariables),
       col=c("black", colors[1:(covariable_num+1)]),
       text.col = c("black", colors[1:(covariable_num+1)]),
       pch = 19, cex = 1.2)

dev.off()

