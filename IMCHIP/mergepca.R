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
nickname=as.character(args[3])
overall_title=as.character(args[4])
covariablefile=as.character(args[5])
colorfile=as.character(args[6])

covariables <- readLines(covariablefile)
colors <- readLines(colorfile)
covariable_num =length(covariables)
#colors <- c("red","blue","green","purple","orange","magenta","yellow")
#covariables <- read.table(file="../covariables.list")
hapmap_color='chartreuse4'

draw_covariable <- function(covariable, covariables, data, colors){
  value <- as.numeric(which(covariables==covariable)[1])
  cv <-as.numeric(data$V9)==value
  points(data[cv,7], data[cv,8], col=colors[value])
}


draw_covariables_ps <- function(covariable, covariables, data, data.hapmap, colors){
  value <- as.numeric(which(covariables==covariable)[1])
  subtable <- data[which(data$V9==value),]
  if (nrow(subtable)>0){
    cov<- gsub(" ","_",covariable)
    cov<- gsub("/","slash",cov)
    filename <- paste0(nickname,"pca_",cov,".ps")
    print(filename)
  postscript(file=filename, paper="letter", horizontal=T)
  
  plot(data.hapmap[,7], data.hapmap[,8], type="p", xlab="PC1", ylab="PC2", main = paste("HapMapIII +",covariable), col = "black")
  u <- par("usr")
  legend(u[1], u[4], xjust=0, yjust=1, 
         legend=c("HapMap",covariable),
         col=c("black", "purple"),
         text.col = c("black", "purple"),
         pch = 19, cex = 1.2)
  
  points(subtable[,7],subtable[,8],col="purple")
  dev.off()
  }
  print("drawn")
}


draw_covariable_sep2 <- function(covariable, covariables, data, data.hapmap, colors){
  plot(data.hapmap[,7], data.hapmap[,8], type="p", xlab="PC1", ylab="PC2", main = paste("HapMapIII +",overall_title,"-",covariable), col = "black")
  value <- as.numeric(which(covariables==covariable)[1])
  subtable <- data[which(data$V9==value),]
  #points(subtable[,7],subtable[,8],col=colors[value])
  points(subtable[,7],subtable[,8],col="purple")
}

draw_covariable_sep <- function(covariable, covariables, data, data.hapmap, colors){
    plot(data.hapmap[,7], data.hapmap[,8], type="p", xlab="PC1", ylab="PC2", main = paste("HapMapIII +",overall_title,"-",covariable), col = "black")
    value <- as.numeric(which(covariables==covariable)[1])
    #subtable=data[which(data$V9==value),]
    #points(subtable[,7],subtable[,8],col=colors[value])
    cv <-as.numeric(data$V9)==value
    points(data[cv,7], data[cv,8], col=colors[value])
  }

#data.hapmap <- read.table(file="hapmappc.txt", header=F, na.string="X")
#data <- read.table(file="projpc.txt", header=F, na.string="X")
data.hapmap <- read.table(file=hapmap_loc, header=F, na.string="X")
pca.data <- read.table(file=table_loc, header=F, na.string="X")

postscript(file=paste0(nickname,"pca.ps"), paper="letter", horizontal=T)
plot(data.hapmap[,7], data.hapmap[,8], type="p", xlab="PC1", ylab="PC2", main = paste("HapMapIII +",overall_title), col = "black")
lapply(covariables, draw_covariable,covariables,pca.data, colors)
u <- par("usr")


legend(u[1], u[4], xjust=0, yjust=1, 
legend=c("HapMap",covariables),
  col=c("black", colors[1:(covariable_num+1)]),
  text.col = c("black", colors[1:(covariable_num+1)]),
  pch = 19, cex = 1.2)

dev.off()

#postscript(file=paste0(nickname,"pca2.ps"), paper="letter", horizontal=T)
#plot(data.hapmap[,7], data.hapmap[,8], type="p", xlab="PC1", ylab="PC2", main = paste("HapMapIII +",overall_title), col = "black")
#points(pca.data[,7],pca.data[,8],col=colors[pca.data[,9]])
#u <- par("usr")


#legend(u[1], u[4], xjust=1, yjust=1, 
#       legend=c("HapMap",covariables),
#       col=c("black", colors[1:(covariable_num+1)]),
#       text.col = c("black", colors[1:(covariable_num+1)]),
#       pch = 19, cex = 1.2)

#dev.off()


postscript(file=paste0(nickname,"pca_covariables.ps"), paper="letter", horizontal=T)
par(mfrow=c(2,2))
lapply(covariables, draw_covariable_sep,covariables,pca.data, data.hapmap, colors)
dev.off()

#postscript(file=paste0(nickname,"pca_covariables2.ps"), paper="letter", horizontal=T)
#par(mfrow=c(2,2))
#lapply(covariables, draw_covariable_sep2,covariables,pca.data, data.hapmap, colors)
#dev.off()


lapply(covariables,draw_covariables_ps,covariables,pca.data,data.hapmap,colors)
