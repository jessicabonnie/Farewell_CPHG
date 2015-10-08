library(ggplot2)
library(dplyr)
print("#########################################################################")
print("#")
print("# GWAS-Pipedream - Standard GWAS QC Pipeline Package")
print("# (c) 2014-2015 JBonnie")
print("# Usage: PCA_PLOTTER_JB.R <pca-projection file> <Study> <chip> <output name for psfile> <xmin> <xmax> <ymin> <ymax> ")
print("#")
print("#########################################################################")

args = commandArgs(TRUE)
print(args)
datafile <- as.character(args[1])
study <- as.character(args[2])
chip <- as.character(args[3])
#covfile <- as.character(args[4])
graphout <- as.numeric(args[4])
boundaries <- FALSE


if (length(args >4)){
  xmin <- as.numeric(args[5])
  xmax <- as.numeric(args[6])
  ymin <- as.numeric(args[7])
  ymax <- as.numeric(args[8])
  boundaries <- TRUE
}

#outfolder <- dirname(datafile)


graphdata <- read.table(file=datafile, header=F, na.string="X")
colnames(graphdata)[7:8] <- c("PC1","PC2")
colnames(graphdata)[colcolumn] <- "Population"
graphdata$Population <- factor(graphdata$Population)


#This counts on Wei-Min's short cut of altering the status before merging with hapmap



#g <- ggplot(mapping=aes(x=PC1,y=PC2)) + geom_point(data=graphdata, aes(group=factor(Population),color = factor(Population))) 

g <- ggplot(mapping=aes(x=PC1,y=PC2)) + geom_point(data=subset(graphdata,V6!=2), aes(group=factor(Population),color = factor(Population)))  + scale_colour_discrete(name  ="HapMapIII")
#g <- g + geom_point(data=subset(graphdata,V6==2),aes(shape=factor(Population)))+ scale_shape_discrete(name="Study")

g <- g + ggtitle(paste('HapMapIII +', study,'-', "\nHumanCoreExome",chip)) 
if (boundaries){
  g <- g + geom_vline(aes(xintercept=xmin),color="#00CC33") + geom_text(aes(x=xmin, label=xmin,y=.02), colour="blue", angle=90, text=element_text(size=11)) 
  g <- g + geom_vline(aes(xintercept=xmax), color="#00CC33") + geom_text(aes(x=xmax, label=xmax,y=.02), colour="blue", angle=90, text=element_text(size=11)) 
  g <- g + geom_hline(aes(yintercept=ymin), color="#00CC33")  + geom_text(aes(x=-.02, label=ymin,y=ymin), colour="blue", text=element_text(size=11)) 
  g <- g+ geom_hline(aes(yintercept=ymax), color="#00CC33") + geom_text(aes(x=-.02, label=ymax,y=ymax), colour="blue", text=element_text(size=11)) 

}

g <- g + geom_point(data=subset(graphdata,V6==2),aes(shape=factor(Population)))+ scale_shape_discrete(name="Study")

postscript(file=paste0(graphout,".ps"), paper="letter", horizontal=T)
g
dev.off()


