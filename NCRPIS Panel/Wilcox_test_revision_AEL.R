library(MASS)
setwd()

#-----Ames traits------
#Look at sweet corn, which is in the Romay folders
dataPath.1="Ames_Data/Conduct_Kchr_Romay/"
dataPath.2="Ames_Data/Pheno+Genotypic_Data_Romay/"

#traitnumber <- c(1,15,16) #1=dT3, 15=dT3/(gT3+aT3), 16=dT3/gT3
traitnames <- c("Sweet")
wilcox <- function(i){
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.",traitnames[i],".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.",traitnames[i],".GWAS.Results.csv",sep=""),head=TRUE)
  
  results.Kchr<-results.Kchr[which((results.Kchr[,2]==4) & (results.Kchr[,3]>=28836485) & (results.Kchr[,3]<=42805440)),]
  results.class<-results.class[which((results.class[,2]==4) & (results.class[,3]>=28836485) & (results.class[,3]<=42805440)),]
  
  combined <- merge(results.Kchr,results.class,by.x="SNP",by.y="SNP")
  wil <- wilcox.test(combined[,4], combined[,14], paired=TRUE) 
  return(wil) 
}

##Error in wilcox.test(combined[, 4], combined[, 14], paired = TRUE) : 
#unused arguments (combined[, 14], paired = TRUE)


result1 <- wilcox(1)

####
#Look at days to silking, which is in the Romay folders
dataPath.1="Ames_Data/Conduct_Kchr_Romay/"
dataPath.2="Ames_Data/Pheno+Genotypic_Data_Romay/"
traitnames <- c("GDD_DTS")

wilcox.test.local <- function(i){
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.",traitnames[i],".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.",traitnames[i],".GWAS.Results.csv",sep=""),head=TRUE)
  
  results.Kchr<-results.Kchr[which((results.Kchr[,2]==8) & (results.Kchr[,3]>=123251085) & (results.Kchr[,3]<=132940001)),]
  results.class<-results.class[which((results.class[,2]==8) & (results.class[,3]>=123251085) & (results.class[,3]<=132940001)),]
  
  combined <- merge(results.Kchr,results.class,by.x="SNP",by.y="SNP")
  wilcox <- wilcox.test(combined[,4],combined[,14], paired=TRUE)  
  return(wilcox) 
}

result2 <- wilcox(1)

#Look at the three traits (days to anthesis, ear height, and plant height) from Peiffer et al. 2014, which is in the Peiffer folders
dataPath.1="Ames_Data/Conduct_Kchr_Peiffer/"
dataPath.2="Ames_Data/Pheno+Geno_Data_Peiffer/"

traitnames <- c("ALL_DTA_BLUE", "ALL_EHT_BLUE", "ALL_PHT_BLUE")
wilcox.test.local <- function(i){
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.",traitnames[i],".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.",traitnames[i],".GWAS.Results.csv",sep=""),head=TRUE)
  
  results.Kchr<-results.Kchr[which((results.Kchr[,2]==8) & (results.Kchr[,3]>=123251085) & (results.Kchr[,3]<=132940001)),]
  results.class<-results.class[which((results.class[,2]==8) & (results.class[,3]>=123251085) & (results.class[,3]<=132940001)),]
  
  merge <- merge(results.Kchr,results.class,by.x="SNP",by.y="SNP")
  wilcox <- wilcox.test(merge[,4], merge[,14], paired=TRUE) 
  return(wilcox) 
}

result3 <- wilcox(1)
result4 <- wilcox(2)
result5 <- wilcox(3)
