print("#########################################################################")
print("#")
print("# GWAS-Pipedream - Standard GWAS QC Pipeline Package")
print("# (c) 2014-2015 JBonnie, WMChen")
print("#")
print("#########################################################################")


args = commandArgs(TRUE)
print(args)
graph_loc=as.character(args[2])
overall_title=gsub("_"," ",as.character(args[3]))
table_loc=as.character(args[1])
covariablefile=as.character(args[5])
nickname=as.character(args[4])
colorfile=as.character(args[6])

covariables <- readLines(covariablefile)
covariable_num =length(covariables)
colors <- readLines(colorfile)
data.raw <- read.table(file = table_loc, header=T)

ysnp_count<-max(data.raw$N_ySNP)
sixthy<-ysnp_count/6
fivesixthy<-5*sixthy
thirdy<-ysnp_count/3
twothirdy<-2*thirdy
half_count <- ysnp_count/2



postscript(graph_loc, paper="letter", horizontal=T)
#jpeg("gender.jpg",width=800, height=800)

N <- sum((data.raw$SEX==1 & data.raw$N_ySNP < half_count) | (data.raw$SEX==2 & data.raw$N_ySNP>half_count) )
#N <- sum((data.raw$SEX==1 & data.raw$N_ySNP < 700) | (data.raw$SEX==2 & data.raw$N_ySNP>700) )
text <- paste(sep="", "Gender Checking in ", overall_title, " Samples (", N, " Samples Mislabeled)")
plot(data.raw$N_ySNP[data.raw$SEX==2], data.raw$xHeterozygosity[data.raw$SEX==2], type="p",
            col = "red", xlim=c(0,max(data.raw$N_ySNP)), ylim=c(0,max(data.raw$xHeterozygosity)),
            main = text,
            xlab="# Y-Chr SNPs", ylab = "X-Chr Heterozygosity")
points(data.raw$N_ySNP[data.raw$SEX==1], data.raw$xHeterozygosity[data.raw$SEX==1], col = "blue")
points(data.raw$N_ySNP[data.raw$SEX==2 & data.raw$N_ySNP>half_count], data.raw$xHeterozygosity[data.raw$SEX==2 & data.raw$N_ySNP>half_count], col = "red")
#points(data.raw$N_ySNP[data.raw$SEX==2 & data.raw$N_ySNP>700], data.raw$xHeterozygosity[data.raw$SEX==2 & data.raw$N_ySNP>700], col = "red")
points(data.raw$N_ySNP[data.raw$SEX==0], data.raw$xHeterozygosity[data.raw$SEX==0], col = "black")

u <- par("usr")
legend(u[2], u[4], xjust=1, yjust=1, 
legend=c("Female", "Male", "Unknown"), col=c("red", "blue", "black"),
text.col = c("red", "blue", "black"), pch = 19, cex = 1.2)

dev.off()


postscript(paste0(nickname,"gender_lines.ps"), paper="letter", horizontal=T)

N <- sum((data.raw$SEX==1 & data.raw$N_ySNP < half_count) | (data.raw$SEX==2 & data.raw$N_ySNP>half_count) )
#N <- sum((data.raw$SEX==1 & data.raw$N_ySNP < 700) | (data.raw$SEX==2 & data.raw$N_ySNP>700) )
text <- paste(sep="", "Gender Checking in ", overall_title, " Samples (", N, " Samples Mislabeled)")
plot(data.raw$N_ySNP[data.raw$SEX==2], data.raw$xHeterozygosity[data.raw$SEX==2], type="p",
     col = "red", xlim=c(0,max(data.raw$N_ySNP)), ylim=c(0,max(data.raw$xHeterozygosity)),
     main = text,
     xlab="# Y-Chr SNPs", ylab = "X-Chr Heterozygosity")
points(data.raw$N_ySNP[data.raw$SEX==1], data.raw$xHeterozygosity[data.raw$SEX==1], col = "blue")
points(data.raw$N_ySNP[data.raw$SEX==2 & data.raw$N_ySNP>half_count], data.raw$xHeterozygosity[data.raw$SEX==2 & data.raw$N_ySNP>half_count], col = "red")
#points(data.raw$N_ySNP[data.raw$SEX==2 & data.raw$N_ySNP>700], data.raw$xHeterozygosity[data.raw$SEX==2 & data.raw$N_ySNP>700], col = "red")
points(data.raw$N_ySNP[data.raw$SEX==0], data.raw$xHeterozygosity[data.raw$SEX==0], col = "black")

abline(v = half_count, col = "green", lty = 3)
abline(v = thirdy, col = "gold", lty = 3)
abline(v = twothirdy, col = "gold", lty = 3)
abline(v = sixthy, col = "violet", lty = 3)
abline(v = sixthy, col = "violet", lty = 3)
abline(v = fivesixthy, col = "violet", lty = 3)
abline(h = 0.05, col = "violet", lty = 3)

u <- par("usr")
legend(u[2], u[4], xjust=1, yjust=1, 
       legend=c("Female", "Male", "Unknown"), col=c("red", "blue", "black"),
       text.col = c("red", "blue", "black"), pch = 19, cex = 1.2)

dev.off()

covariable_gender_sep <-function(covariable,covariables,data){
  value=as.numeric(which(covariables==covariable)[1])
  subtable <- data[which(as.numeric(data$Covariable)==value),]
  if (nrow(subtable) > 2){
    cov_as_filename <- gsub("/","slash",gsub(" ","_",covariable))
    postscript(paste0(nickname,"gender_",cov_as_filename,".ps"), paper="letter", horizontal=T)
    subysnp_count<-max(subtable$N_ySNP)
    subhalf_count <- subysnp_count/2
    subN = sum((subtable$SEX==1 & subtable$N_ySNP < subhalf_count) | (subtable$SEX==2 & subtable$N_ySNP>subhalf_count) )
    subtext <- paste(sep="", "Gender Checking in ", covariable, " Samples (", subN, " Samples Mislabeled)")
    plot(subtable$N_ySNP[subtable$SEX==2], subtable$xHeterozygosity[subtable$SEX==2], type="p",
         col = "red", xlim=c(0,max(data$N_ySNP)), ylim=c(0,max(data$xHeterozygosity)),
         main = subtext,
         xlab="# Y-Chr SNPs", ylab = "X-Chr Heterozygosity")
    points(subtable$N_ySNP[subtable$SEX==1], subtable$xHeterozygosity[subtable$SEX==1], col = "blue")
    points(subtable$N_ySNP[subtable$SEX==2 & subtable$N_ySNP>subhalf_count], subtable$xHeterozygosity[subtable$SEX==2 & subtable$N_ySNP>subhalf_count], col = "red")
    #points(data.raw$N_ySNP[data.raw$SEX==2 & data.raw$N_ySNP>700], data.raw$xHeterozygosity[data.raw$SEX==2 & data.raw$N_ySNP>700], col = "red")
    points(subtable$N_ySNP[subtable$SEX==0], subtable$xHeterozygosity[subtable$SEX==0], col = "black")
  
  
    abline(v = half_count, col = "green", lty = 3)
    abline(v = thirdy, col = "gold", lty = 3)
    abline(v = twothirdy, col = "gold", lty = 3)
    abline(v = sixthy, col = "violet", lty = 3)
    abline(v = sixthy, col = "violet", lty = 3)
    abline(v = fivesixthy, col = "violet", lty = 3)
    abline(h = 0.05, col = "violet", lty = 3)
  
    u <- par("usr")
    legend(u[2], u[4], xjust=1, yjust=1, 
           legend=c("Female", "Male", "Unknown"), col=c("red", "blue", "black"),
           text.col = c("red", "blue", "black"), pch = 19, cex = 1.2)
    
    dev.off()
  }
}

covariable_colored <-function(covariable,covariables,colors,data){
  value=as.numeric(which(covariables==covariable)[1])
  subtable <- data[which(as.numeric(data$Covariable)==value),]
  points(subtable$N_ySNP,subtable$xHeterozygosity,col=colors[value])
    
  
}

lapply(covariables,covariable_gender_sep,covariables,data.raw)

postscript(paste0(nickname,"gender_covariable.ps"), paper="letter", horizontal=T)
plot(data.raw$N_ySNP, data.raw$xHeterozygosity, type="p",
     col = "transparent", xlim=c(0,max(data.raw$N_ySNP)), ylim=c(0,max(data.raw$xHeterozygosity)),
     main = paste(sep="", "Gender in ", overall_title, " Samples, Colored by Covariable"),
     xlab="# Y-Chr SNPs", ylab = "X-Chr Heterozygosity")

abline(v = half_count, col = "green", lty = 3)
abline(v = thirdy, col = "gold", lty = 3)
abline(v = twothirdy, col = "gold", lty = 3)
abline(v = sixthy, col = "violet", lty = 3)
abline(v = sixthy, col = "violet", lty = 3)
abline(v = fivesixthy, col = "violet", lty = 3)
abline(h = 0.05, col = "violet", lty = 3)

legend(u[2], u[4], xjust=1, yjust=1, 
       legend=covariables, col=colors[1:(covariable_num+1)],
       text.col = colors[1:(covariable_num+1)], pch = 19, cex = 1.2)
lapply(covariables,covariable_colored,covariables,colors,data.raw)

dev.off()
