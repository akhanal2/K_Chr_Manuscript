setwd()
source("http://peterhaschke.com/Code/multiplot.R")
#***************************************************************************
#                   Specific traits & specific genomic regions
#***************************************************************************

multiplot(p1, p2, p3, p4, cols=2)

library(gridExtra)
library(ggplot2)


#-----Ames traits------
#Look at sweet corn, which is in the Romay folders
dataPath.1="Conduct_Kchr_Romay/"
dataPath.2="Pheno+Genotypic_Data_Romay/"

#traitnumber <- c(1,15,16) #1=dT3, 15=dT3/(gT3+aT3), 16=dT3/gT3
traitnames <- "Sweet"
boxplot <- function(i){
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.",traitnames[i],".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.",traitnames[i],".GWAS.Results.csv",sep=""),head=TRUE)
 
  results.Kchr<-results.Kchr[which((results.Kchr[,2]==4) & (results.Kchr[,3]>=28836485) & (results.Kchr[,3]<=42805440)),]
  results.class<-results.class[which((results.class[,2]==4) & (results.class[,3]>=28836485) & (results.class[,3]<=42805440)),]
  
  results.Kchr$Method<-rep("K_chr", times = nrow(results.Kchr)) 
  results.class$Method<-rep("Trad MLM", times = nrow(results.class)) 
  results.Kchr <-results.Kchr[,c(4,ncol(results.Kchr))] #use raw p value
  results.class <-results.class[,c(4,ncol(results.class))]
  combined <- rbind(results.Kchr,results.class)
  colnames(combined)<-c("p_value","Method")
  
  tmp.plot <-ggplot(combined,
                    aes(x=Method, y=-log10(p_value), fill=Method)) +
    geom_boxplot() + 
    guides(fill=FALSE) + 
    stat_summary(fun.y=mean, geom="point", shape=5, size=2) 
  # + ggtitle(paste("ZmVTE1_", traitnames[i],sep=""))
  return(tmp.plot) 
}

boxplot(1) 

####
#Look at days to silking, which is in the Romay folders
dataPath.1="Conduct_Kchr_Romay/"
dataPath.2="Pheno+Genotypic_Data_Romay/"

#traitnumber <- c(1,15,16) #1=dT3, 15=dT3/(gT3+aT3), 16=dT3/gT3
traitnames <- "GDD_DTS"
boxplot <- function(i){
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.",traitnames[i],".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.",traitnames[i],".GWAS.Results.csv",sep=""),head=TRUE)
  
  results.Kchr<-results.Kchr[which((results.Kchr[,2]==8) & (results.Kchr[,3]>=123251085) & (results.Kchr[,3]<=132940001)),]
  results.class<-results.class[which((results.class[,2]==8) & (results.class[,3]>=123251085) & (results.class[,3]<=132940001)),]
  
  results.Kchr$Method<-rep("K_chr", times = nrow(results.Kchr)) 
  results.class$Method<-rep("Trad MLM", times = nrow(results.class)) 
  results.Kchr <-results.Kchr[,c(4,ncol(results.Kchr))] #use raw p value
  results.class <-results.class[,c(4,ncol(results.class))]
  combined <- rbind(results.Kchr,results.class)
  colnames(combined)<-c("p_value","Method")
  
  tmp.plot <-ggplot(combined,
                    aes(x=Method, y=-log10(p_value), fill=Method)) +
    geom_boxplot() + 
    guides(fill=FALSE) + 
    stat_summary(fun.y=mean, geom="point", shape=5, size=2) 
  # + ggtitle(paste("ZmVTE1_", traitnames[i],sep=""))
  return(tmp.plot) 
}

boxplot(1) 

#Look at the three traits from Peiffer et al. 2014, which is in the Peiffer folders
dataPath.1="Conduct_Kchr_Peiffer/"
dataPath.2="Pheno+Geno_Data_Peiffer/"


traitnames <- c("ALL_DTA_BLUE", "ALL_EHT_BLUE", "ALL_PHT_BLUE")
boxplot(1) 
boxplot(2)
boxplot(3)




