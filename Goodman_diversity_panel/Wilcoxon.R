library(MASS)


#-----Tocochromanol-Specific Traits-ZmVTE1------
dataPath.1="K_chr/Tocochromanol/Combined_Tocochromanol_Kchr/"
dataPath.2="Classical/Tocochromanol/Combined_Tocochromanol_classical/"

traitnumber <- c(1,15) #1=dT3, 15=dT3/(gT3+aT3)
traitnames <- c("dT3", "dT3/(gT3+aT3)")
wilcox <- function(i){
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.traittrans",traitnumber[i],".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.traittrans",traitnumber[i],".GWAS.Results.csv",sep=""),head=TRUE)
  
  results.Kchr<-results.Kchr[which((results.Kchr[,2]==5) & (results.Kchr[,3]>=132172312) & (results.Kchr[,3]<=136200852)),]
  results.class<-results.class[which((results.class[,2]==5) & (results.class[,3]>=132172312) & (results.class[,3]<=136200852)),]
  
  combined <- merge(results.Kchr,results.class,by.x="SNP",by.y="SNP")
  wil <- wilcox.test(combined[,4], combined[,14], paired=TRUE) 
  return(wil) 
}

##Error in wilcox.test(combined[, 4], combined[, 14], paired = TRUE) : 
#unused arguments (combined[, 14], paired = TRUE)

wilcox(1) 
wilcox(2)

#-----Tocochromanol-Specific Traits-ZmVTE4------
dataPath.1="Combined_Tocochromanol_Kchr/"
dataPath.2="Combined_Tocochromanol_classical/"

traitnumber <- c(6,13,14,19) #6=aT, 13=dT/aT, 14=gT/(gT+aT), 19=aT/gT

traitnames <- c("aT", "dT/aT", "gT/(gT+aT)", "aT/gT")
wilcox.test <- function(i){
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.traittrans",traitnumber[i],".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.traittrans",traitnumber[i],".GWAS.Results.csv",sep=""),head=TRUE)
  
  results.Kchr<-results.Kchr[which((results.Kchr[,2]==5) & (results.Kchr[,3]>=197512964) & (results.Kchr[,3]<=204549958)),]
  results.class<-results.class[which((results.class[,2]==5) & (results.class[,3]>=197512964) & (results.class[,3]<=204549958)),]
  
  merge <- merge(results.Kchr,results.class,by.x="SNP",by.y="SNP")
  wilcox <- wilcox.test(merge[,4], merge[,14], paired=TRUE)  
  return(wilcox) 
}

wilcox.test(1) 
wilcox.test(2)
wilcox.test(3) 
wilcox.test(4)

#-----Carotenoids-Specific Traits-ZEP1------
dataPath.1="Combined_Carotenoids_Kchr/"
dataPath.2="Combined_Carotenoids_classical/"

traitnumber <- c(2,11,14) #2=Zeaxanthin, 11=TbetaXantho, 14=b2aXantho
traitnames <- c("Zeaxanthin", "TbetaXantho", "b2aXantho")
wilcox.test <- function(i){
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.traittrans",traitnumber[i],".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.traittrans",traitnumber[i],".GWAS.Results.csv",sep=""),head=TRUE)
  
  results.Kchr<-results.Kchr[which((results.Kchr[,2]==2) & (results.Kchr[,3]>=44035247) & (results.Kchr[,3]<=45199441)),]
  results.class<-results.class[which((results.class[,2]==2) & (results.class[,3]>=44035247) & (results.class[,3]<=45199441)),]
  
  merge <- merge(results.Kchr,results.class,by.x="SNP",by.y="SNP")
  wilcox <- wilcox.test(merge[,4], merge[,14], paired=TRUE) 
  return(wilcox) 
}

wilcox.test(1) 
#wilcox.test(2) trait11 can't be run because no original file exist 
wilcox.test(3)


#-----Carotenoids-Specific Traits-LUT1------
#library(gridExtra)
#library(ggplot2)
dataPath.1="Combined_Carotenoids_Kchr/"
dataPath.2="Combined_Carotenoids_classical/"

traitnumber <- c(3,17,30) #3=Zeinoxanthin, 17=A2zeino, 30=zeino2lut
traitnames <- c("Zeinoxanthin", "A2zeino", "zeino2lut")
wilcox.test <- function(i){
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.traittrans",traitnumber[i],".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.traittrans",traitnumber[i],".GWAS.Results.csv",sep=""),head=TRUE)
  
  results.Kchr<-results.Kchr[which((results.Kchr[,2]==1) & (results.Kchr[,3]>=86586169) & (results.Kchr[,3]<=87600028)),]
  results.class<-results.class[which((results.class[,2]==1) & (results.class[,3]>=86586169) & (results.class[,3]<=87600028)),]
  
  merge <- merge(results.Kchr,results.class,by.x="SNP",by.y="SNP")
  wilcox <- wilcox.test(merge[,4], merge[,14], paired=TRUE) 
  return(wilcox) 
}

wilcox.test(1)
wilcox.test(2)
wilcox.test(3)