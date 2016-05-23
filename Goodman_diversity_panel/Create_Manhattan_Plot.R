# GWAS with 4K SNPs and 55K SNPs   
#Set the working directory
setwd("/Users/adminuser/Box Sync/Data_for_Angela/Manuscript/Create_Manhattan_Plot")
home.dir <- getwd()

#Source in the appropriate files
library('MASS')
source("http://zzlab.net/GAPIT/emma.txt")
source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")


#Select a subdirectory
sub.dir <- "Carotenoid_Traits_in_Novel_Genomic_Regions"
setwd(paste(home.dir, "/", sub.dir, sep = ""))
trait.class <- "Carotenoid"

#Read in the data
the.data <- read.table("Input_File_Novel_Manhattan_Plot_Carotenoids.txt", head = TRUE)
#Get rid of the first column so that this can be read properly into GAPIT.Manhattan
the.data <- the.data[,-1]

#Make the Manhattan plot
GAPIT.Manhattan(GI.MP = the.data, name.of.trait = trait.class)







##############################################################################################
GAPIT.Manhattan <- function(GI.MP = NULL, name.of.trait = "Trait", 
                            plot.type = "Genomewise",DPP=50000,cutOff=0.01){
  #Object: Make a Manhattan Plot
  #Options for plot.type = "Separate_Graph_for_Each_Chromosome" and "Same_Graph_for_Each_Chromosome" 
  #Output: A pdf of the Manhattan Plot
  #Authors: Alex Lipka, Zhiwu Zhang, and Meng Li 
  # Last update: May 10, 2011 
  ##############################################################################################
  
  print("Manhattan ploting...")
  
  #do nothing if null input
  if(is.null(GI.MP)) return
  
  GI.MP=matrix(as.numeric(as.matrix(GI.MP) ) ,nrow(GI.MP),ncol(GI.MP))
  
  #Remove all SNPs that do not have a choromosome and bp position
  GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
  GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
  
  #Remove all SNPs that have P values above 0 (not na etc)
  #GI.MP <- GI.MP[GI.MP[,3]>0,]
  numMarker=nrow(GI.MP)
  #bonferroniCutOff=-log10(cutOff/numMarker)
  
  print("debug for plot")
  print(numMarker)
  #print(bonferroniCutOff)
  
  #Replace P the -log10 of the P-values
  y.lim.vector <- NULL
  for(i in 3:ncol(GI.MP)){
    GI.MP[,i] <-  -log10(GI.MP[,i])
    y.lim.vector <- c(y.lim.vector, ceiling(max(GI.MP[which(!is.na(GI.MP[,i])),i])))
   }#end for(i in 3:ncol(GI.MP))
  
  #set the limit of the Y-axis
  y.lim <- max(y.lim.vector)
  
  
  chm.to.analyze <- unique(GI.MP[,1]) 
  chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
  numCHR= length(chm.to.analyze)
  
 
  
  #Genomewise plot
  print("Manhattan ploting Genomewise")
  
  GI.MP <- GI.MP[order(GI.MP[,2]),]
  GI.MP <- GI.MP[order(GI.MP[,1]),]
  #color.vector <- rep(c("orangered","navyblue"),numCHR)
  ticks=NULL
  lastbase=0
  
  #print("Manhattan data sorted")
  #print(chm.to.analyze) 
  
  #change base position to accumulatives
  for (i in chm.to.analyze)
  {
    index=(GI.MP[,1]==i)
    ticks <- c(ticks, lastbase+mean(GI.MP[index,2])) 
    GI.MP[index,2]=GI.MP[index,2]+lastbase
    lastbase=max(GI.MP[index,2])
  }
  
  ticks.end <- NULL
  for(j in chm.to.analyze[-length(chm.to.analyze)]){
    data.subset <- GI.MP[which(GI.MP[,1]==j),]
    ticks.end <- c(ticks.end, max(data.subset[,2]))
  }
  
  #print("Manhattan chr processed")
  
  x0 <- as.numeric(GI.MP[,2])
  y0 <- as.numeric(GI.MP[,3])
  z0 <- as.numeric(GI.MP[,1])
  position=order(y0,decreasing = TRUE)
  index0=GAPIT.Pruning(y0[position],DPP=DPP)
  index=position[index0]
  x=x0
  y=y0
  z=z0
  
  #print("Manhattan XY created")
  
  #Pick a vector of reasonable colors of size 10. If this code ever needs to be
  # generalized, then this section of the code (creation of a vector of
  # reasonable colors) will need to be revised.
  color.vector <- c("orange", "red", "purple", "black", "navyblue", "gray42",
                    "orangered1", "violetred1", "azure4", "dodgerblue")
  
  pdf(paste("Novel.Genomic.Regions.", name.of.trait,".Manhattan-Plot.pdf" ,sep = ""), width = 10)
  par(mar = c(7,7,7,1))
  
  plot(y~x,xlab=expression(Chromosome),ylab=expression(-log[10](italic(p))) ,
       cex.lab=2,col="blue",axes=FALSE,type = "p",pch=20, ylim = c(0, y.lim))
  #abline(h=bonferroniCutOff,col="forestgreen")
  if(ncol(GI.MP) >= 4) {
    for(i in 4:ncol(GI.MP)){
      y <- as.numeric(GI.MP[,i])
     # par(new = TRUE)
      lines(y~x,type = "p",pch=20, col = color.vector[i-3], ylim = c(0, y.lim))
    }#end for(i in 4:ncol(GI.MP))
  }#end if(ncol(GI.MP) > 2)
  
  
  
  
  axis(1, at=ticks,cex.axis=1.5,labels=chm.to.analyze,tick=F)
  axis(2, at=1:y.lim,cex.axis=1.5,labels=1:y.lim,tick=F)
  abline(v = ticks.end, col = "lightgrey", lwd = 1.5, cex = 1.0)
  box()
  dev.off()
  print("Manhattan done Genomewise")
  


  print("GAPIT.Manhattan accomplished successfully!")
} #end of GAPIT.Manhattan










