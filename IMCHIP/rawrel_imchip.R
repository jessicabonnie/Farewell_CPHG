print("#########################################################################")
print("#")
print("# GWAS-Pipedream - Standard GWAS QC Pipeline Package")
print("# (c) 2014-2015 JBonnie, WMChen")
print("#")
print("#########################################################################")

args = commandArgs(TRUE)
print(args)
table_loc=as.character(args[1])
crosspair_loc=as.character(args[2])
nickname=as.character(args[3])
overall_title=gsub("_"," ",as.character(args[4]))
covariablefile=as.character(args[5])
colorfile=as.character(args[6])
covariables <- readLines(covariablefile)
covariable_num =length(covariables)
colors <- readLines(colorfile)

kinship_min=0.02
SDthresh=0.1
FDthresh=.18

draw_covariable_pairs <- function(covariable, covariables, data, colors){
  value=as.numeric(which(covariables==covariable)[1])
  c1v <-as.numeric(data$Covariable1)==value & as.numeric(data$Kinship) > kinship_min
  c2v <-as.numeric(data$Covariable2)==value & as.numeric(data$Kinship) > kinship_min
  points(data[c1v,]$IBS0, data[c1v,]$Kinship, col=colors[value])
  points(data[c2v,]$IBS0, data[c2v,]$Kinship, col=colors[value],pch=4)
}
draw_intra_covariables <- function(covariable, covariables, data, colors){
  value=as.numeric(which(covariables==covariable)[1])
  cv <-as.numeric(data$Covariable1)==value & as.numeric(data$Covariable2)==value & as.numeric(data$Kinship) > kinship_min
  points(data[cv,]$IBS0, data[cv,]$Kinship, col=colors[value])
}

draw_covariable_ps <- function(covariable,covariables,data,colors){
  cov_as_filename <- gsub("/","slash",gsub(" ","_",covariable))
  postscript(paste0(nickname,"rawrel_",cov_as_filename,".ps"), paper="letter", horizontal=T)
  value <- as.numeric(which(covariables==covariable)[1])
  subtable <- data[which(as.numeric(data$Covariable1)==value & as.numeric(data$Covariable2)==value & as.numeric(data$Kinship)>kinship_min),]
  print(head(subtable))
  subN <- sum(subtable$Kinship > SDthresh)
  #c1v <-as.numeric(data$Covariable1)==value & as.numeric(data$Covariable2)==value & as.numeric(data$Kinship) > 0.02
  #c2v <-as.numeric(data$Covariable2)==value & as.numeric(data$Kinship) > 0.02
  #plot(subtable$IBS0[subtable$Kinship > 0.02], subtable$Kinship[subtable$Kinship > 0.02], type="p",
  plot(subtable$IBS0, subtable$Kinship, type="p",
       col = "blue", cex.lab=1.3,xlim=c(0,max(data$IBS0)), ylim=c(0,0.5),
       # ylim=c(kinship_min,max(data$Kinship)),
       main = paste("Relationship in", covariable,"Among Individuals (",subN,"pairs with Kinship > ",SDthresh," )"),
       xlab="Proportion of Zero IBS", ylab = "Estimated Kinship Coefficient")
  abline(h = 0.3536, col = "red", lty = 4)
  abline(h = FDthresh, col = "green", lty = 4)
  abline(h = SDthresh, col = "gold", lty = 4)
  

  dev.off()
  #points(data[c2v,]$IBS0, data[c2v,]$Kinship, col=colors[value],pch=4)
}


data <- read.table(file = table_loc, header=T)
crosspair.data <- read.table(file=crosspair_loc,header=T)
N <- sum(data$Kinship > SDthresh)
lapply(covariables, draw_covariable_ps,covariables,data, colors)
postscript(paste0(nickname,"rawrel_covariable_pairs.ps"), paper="letter", horizontal=T)
#jpeg("rawrel_covariable_pairs.jpg",width=800, height=800)
#data <- read.table(file = "dn6.kin0", header=T)

plot(data$IBS0[data$Kinship > kinship_min], data$Kinship[data$Kinship > kinship_min], type="p",
     col = "lightgray", cex.lab=1.3,
     main = paste("Relationship in", overall_title,"Among Individuals (",N,"pairs with Kinship >",SDthresh,")"),
     xlab="Proportion of Zero IBS", ylab = "Estimated Kinship Coefficient", xlim=c(0,max(data$IBS0)),ylim=c(0,0.5))

lapply(covariables, draw_covariable_pairs,covariables,crosspair.data, colors)

u <- par("usr")
legend(u[2], u[4], xjust=1, yjust=1, 
       legend=covariables,
       col=colors[1:(covariable_num+1)],
       text.col = colors[1:(covariable_num+1)],
       pch = 19, cex = 1.2)

abline(h = 0.3536, col = "red", lty = 4)
abline(h = FDthresh, col = "green", lty = 4)
abline(h = SDthresh, col = "gold", lty = 4)


dev.off()


postscript(paste0(nickname,"rawrel.ps"), paper="letter", horizontal=T)
#jpeg("rawrel.jpg",width=800, height=800)
plot(data$IBS0[data$Kinship > kinship_min], data$Kinship[data$Kinship > kinship_min], type="p",
     col = "blue", cex.lab=1.3,
     main = paste("Relationship in", overall_title,"Among Individuals (",N,"pairs with Kinship > 0.1768 )"),
     xlab="Proportion of Zero IBS", ylab = "Estimated Kinship Coefficient",ylim=c(0,0.5),xlim=c(0,max(data$IBS0)))


abline(h = 0.3536, col = "red", lty = 4)
abline(h = FDthresh, col = "green", lty = 4)
abline(h = SDthresh, col = "gold", lty = 4)


dev.off()


postscript(paste0(nickname,"rawrel_covariables.ps"), paper="letter", horizontal=T)
#jpeg("rawrel_covariables.jpg",width=800, height=800)
#data <- read.table(file = "dn6.kin0", header=T)

plot(data$IBS0[data$Kinship > kinship_min], data$Kinship[data$Kinship > kinship_min], type="p",
     col = "lightgray", cex.lab=1.3,
     main = paste("Relationship in", overall_title,"Among Individuals (",N,"pairs with Kinship >",SDthresh,")"),
     xlab="Proportion of Zero IBS", ylab = "Estimated Kinship Coefficient",ylim=c(0,0.5),xlim=c(0,max(data$IBS0)))
lapply(covariables, draw_intra_covariables, covariables,data,colors)
lapply(covariables, draw_covariable_pairs,covariables,crosspair.data, colors)

u <- par("usr")
legend(u[2], u[4], xjust=1, yjust=1, 
       legend=covariables,
       col=colors[1:(covariable_num+1)],
       text.col = colors[1:(covariable_num+1)],
       pch = 19, cex = 1.2)

abline(h = 0.3536, col = "red", lty = 4)
abline(h = FDthresh, col = "green", lty = 4)
abline(h = SDthresh, col = "gold", lty = 4)

dev.off()


lapply(covariables, draw_covariable_ps,covariables,data, colors)