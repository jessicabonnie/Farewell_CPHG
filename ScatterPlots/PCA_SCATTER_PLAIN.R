library(ggplot2)
library(dplyr)
library(reshape2)   # for melt(...)
library(plyr)       # for .(...)
library(data.table)
print("#########################################################################")
print("#")
print("# GWAS-Pipedream - Standard GWAS QC Pipeline Package")
print("# (c) 2014-2015 JBonnie")
print("# It is currently expected that the 27th column will be the one to color by.")
print("# It is currently expected that the 6th column indicates which sample should be drawn in color.
      Samples with a 2 in the status column will be black.")
print("# Usage: PCA_SCATTER_JB.R <pca-projection file> <Study> <chip>  <output file name>  ")
print("#")
print("#########################################################################")

args = commandArgs(TRUE)
print(args)
datafile <- as.character(args[1])
study <- as.character(args[2])
chip <- as.character(args[3])
outfilename <- as.character(args[4])
overall_title="JDRFDN"
#outfolder <- dirname(datafile)


graphdata <- read.table(file=datafile, header=F, na.string="X")
colnames(graphdata) <- c("FID","IID","Father","Mother","Sex","Status",paste0("PC",rep(1:20)),"Population")

graphdata$Population <- factor(graphdata$Population)
#blackname <- unique(as.character(subset(graphdata,Status==2)$Population))
#graphdata <- graphdata[,c(1:16,27)]

#This counts on Wei-Min's short cut of altering the status before merging with hapmap

# study_data$Population <- factor(as.character(study_data$Population))
# 
# 
# graph_pc <- function (col1, col2){
#   pcplot <- ggplot(mapping=aes_string(x=col1,y=col2)) + geom_point(data=subset(graphdata,Status!=2),
#                                                                    aes(group=factor(Population),
#                                                                        color = factor(Population)))  + scale_colour_discrete(name  ="HapMapIII")
#   pcplot <- pcplot + geom_point(data=subset(graphdata,Status==2),aes(shape=factor(Population)))+ scale_shape_discrete(name="Study")
#   return (pcplot)
# }
# 
# 
# 
# library(GGally)
# 
# pcs <-c("PC1","PC2","PC3","PC4")
# pc_scatter <- ggpairs(subset(graphdata,Status!=2)[,c(pcs,"Status","Population")] , columns=1:4,title = paste(blackname, "Principal Components"),
#                       diag="blank",color="Population",legends=TRUE,aes(group=factor(Population),color = factor(Population))) 
# 
# for (i in rep(1:4)){
#   
#   for (j in rep(2:4)){
# 
#     pcx <- paste0("PC",i)
#     pcy <- paste0("PC",j)
#     print(pcx)
#     print(pcy)
#     cp <- graph_pc(pcx,pcy)
#     pc_scatter <-putPlot(pc_scatter, cp, i, j)
#     
#   }
# }
# pc_scatter
# # pc_scatter
# 
# custom_plot <- ggpairs(graphdata[,pcs], upper = "blank", lower="blank", title = paste(blackname, "Principal Components"))
# 


#dev.off()

#pcs <- paste0("PC",rep(1:10))
#pcs <-c("PC1","PC2","PC3","PC4")

graph_these <- function(pclist){
#take only what you need from the table and then add row id values
graph_df <- with(graphdata, data.table(id=1:nrow(graphdata), group=Population, graphdata[,pclist]))
# melt the table (as in create a cell for every single variable and then fill in)
graphtable_melt <- melt(graph_df,id=1:2, variable.name="pcnum1", value.name="pcval1")
setkey(graphtable_melt,id,group)

#make a copy of the melt but give the columns different names, so that we can have PCA values graphed against each other on each axis.
gmelty <- graphtable_melt[,list(pcnum2=pcnum1,pcval2=pcval1),key="id,group"]

# put the two together
combmelt <- graphtable_melt[gmelty,allow.cartesian=T]
setkey(combmelt,pcnum1,pcnum2,group)


# add ranges maxes and mins
combmelt <- combmelt[,list(id, group, pcval1, pcval2, min.x=min(pcval1),min.y=min(pcval2),
                           range.x=diff(range(pcval1)),range.y=diff(range(pcval2))),by="pcnum1,pcnum2"]


#create table when pc1=pc2 for density plot
#d  <- zz[pcnum1==pcnum2,list(x=density(pcval1)$x,
#                   y=min.y+range.y*density(pcval1)$y/max(density(pcval1)$y)),
#         by="pcnum1,pcnum2,group"]

# We only want to graph each PC combination once, so we now filter the table
#combmelt<-filter(combmelt, gsub("PC","",pcnum1) > gsub("PC","",pcnum2))


facetplot <- ggplot(data=combmelt,aes(x=pcval1,y=pcval2))


  facetplot <- facetplot +
    geom_point(aes(x=pcval1, y=pcval2, color=factor(group))) +
    facet_grid(pcnum2~pcnum1, scales="free") +
    ggtitle(paste(overall_title,'--', study,'\n', '\nHumanCoreExome',chip)) + 
    scale_color_discrete(name="") + labs(x="", y="")

facetplot

}

postscript(file=paste0(outfilename,".ps"), paper="letter", horizontal=T)
pclist1  <- paste0("PC",rep(1:5))
graph_these(pclist1)

pclist2  <- paste0("PC",rep(6:10))
graph_these(pclist2)

dev.off()

