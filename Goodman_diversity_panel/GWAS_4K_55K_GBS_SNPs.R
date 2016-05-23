# GWAS with 4K SNPs and 55K SNPs   

#######################################################################################
#Step 1: Import files
#######################################################################################
setwd() 
library('MASS')
source("http://zzlab.net/GAPIT/emma.txt")
source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")
source("mt.rawp2adjp.r")
mydata  <- read.delim("4K_SNPsmdp_genotype_test1.hmp.txt", head = FALSE) #read in genotype 
myY  <- read.table("mdp_traits.txt", head = TRUE) #read in phenotype 

#######################################################################################
#Step 2: Set result directory and run GAPIT
#######################################################################################
for (i in 1:10){
  #looping through a total of 10 chromosomes 
  
  obs <- mydata[which(as.character(mydata[,3])==as.character(i)),] 
  #parse genotype data (extract SNPs from chromosome i)
  
  first.row <- mydata[1,] 
  #grab SNP names  
  
  sub <- rbind(first.row,obs)
  #create a file that has the parse data and corresponding SNP names
    
  subdataPath.1=paste("Kchr_chrom_", as.character(i), sep="")
  setwd(subdataPath.1)
  #set working directory to grab kinship and principal components files
  
  myKI <- read.csv("GAPIT.Kin.Loiselle.csv", head = FALSE)
  #read in kinship 
  
  myCV <- read.csv("GAPIT.PCA.csv", head = TRUE)
  #read in principal components 
  
  subdataPath.2=paste("Kchr_GWAS_4K_chrom_", as.character(i), sep="")
  dir.create(subdataPath.2, showWarnings = FALSE)
  setwd(subdataPath.2)
  #create a new directory to store GWAS results 
  
  myGAPIT <- GAPIT(
    Y=myY,        #phenotype data
    G=sub,      		#genotype data
    file.path=NULL,		#Set the location of genotype files to NULL
    KI=myKI,				#kinship data
    CV=myCV,				#covariate variables of fixed effects, such as population structure
    group.from=272,		#Lower bound for number of group
    group.to=272,			#Upper bound for number of group
    group.by=1,				#rang between 1 and number of individuals, smaller the finner optimization 
    kinship.cluster="average", 			#clustering method: "average", "complete", "ward", "single", "mcquitty", "median","centroid". Example: CA=c("complete","average")
    kinship.group="Mean",     		#Group kinship tppe:  "Mean", "Max", "Min", "Median". Example: KT=c("Mean","Max")
    SNP.P3D=TRUE,		#This is the option to use P3D (TRUE) or not (FALSE)    
    SNP.impute = "Major",
    SNP.MAF = 0.05,
    cutOff = 0.00,
    Model.selection = TRUE   
  )
 }


# GWAS with GBS SNPs
#######################################################################################
#Step 1: Import files
#######################################################################################
setwd() 
library('MASS')
source("http://zzlab.net/GAPIT/emma.txt")
source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")
source("mt.rawp2adjp.r")
myY  <- read.table("mdp_traits.txt", head = TRUE) #read in phenotype  
mydataPath=""

#######################################################################################
#Step 2: Set result directory and run GAPIT
#######################################################################################
for (i in 1:10){
  #looping through a total of 10 chromosomes 
  
  obs <- mydata[which(as.character(mydata[,3])==as.character(i)),] 
  #parse genotype data (extract SNPs from chromosome i)
  
  first.row <- mydata[1,] 
  #grab SNP names  
  
  sub <- rbind(first.row,obs)
  #create a file that has the parse data and corresponding SNP names
  
  subdataPath.1=paste("Kchr_chrom_", as.character(i), sep="")
  setwd(subdataPath.1)
  #set working directory to grab kinship and principal components files
  
  myKI <- read.csv("GAPIT.Kin.Loiselle.csv", head = FALSE)
  #read in kinship 
  
  myCV <- read.csv("GAPIT.PCA.csv", head = TRUE)
  #read in principal components 
  
  
  subdataPath.2=paste("Kchr_GWAS_GBS_chrom_", as.character(i), sep="")
  dir.create(subdataPath.2, showWarnings = FALSE)
  setwd(subdataPath.2)
  #create a new directory to store GWAS results 
  
  myGAPIT <- GAPIT(
    Y=myY,        #phenotype data
    G=NULL,       #genotype data,set it to NULL with multiple genotype files 
    file.path=mydataPath,  	#The location of genotype files
    file.from = i, 
    file.to = i,
    file.total = 1,
    file.G="282_20111217_scv10mF8maf002_mgs_E1pLD5.imp95_1024.c", #Common file name (before the numerical part), set it to NULL if for single file
    file.Ext.G="hmp.txt",
    
    KI=myKI,				#kinship data
    CV=myCV,				#covariate variables of fixed effects, such as population structure
    group.from=272,		#Lower bound for number of group
    group.to=272,			#Upper bound for number of group
    group.by=1,				#rang between 1 and number of individuals, smaller the finner optimization 
    kinship.cluster="average", 			#clustering method: "average", "complete", "ward", "single", "mcquitty", "median","centroid". Example: CA=c("complete","average")
    kinship.group="Mean",     		#Group kinship tppe:  "Mean", "Max", "Min", "Median". Example: KT=c("Mean","Max")
    SNP.P3D=TRUE,		#This is the option to use P3D (TRUE) or not (FALSE)   
    SNP.impute = "Major",
    SNP.MAF = 0.05,
    cutOff = 0.00,
    Model.selection = TRUE   
  )
}

