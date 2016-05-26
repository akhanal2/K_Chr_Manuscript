#Step 1: Set data directory and import files
#######################################################################################
setwd()
home.dir <- getwd()

library('MASS')
library(multtest)
library(gplots)

source("http://zzlab.net/GAPIT/emma.txt")
source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")


#chrom.1 <- mydata[,3]==1
#chrom.10 <- mydata[,3]==10
#sum(chrom.1)
#sum(chrom.10)

myY  <- read.table(paste(home.dir, "/Pheno+Geno_Data/Romay_etal_2013_GenomeBiol_phenotypes-130503_for_R.txt", sep = ""), head = TRUE)
#head(myY)

#Get rid of kernel color because we are not going to analyze this trait
myY <- myY[,-2]

#Create unique phenotype files for each of the two traits
myY.sweet <- myY[,-3]
myY.sweet <- myY.sweet[-which(is.na(myY.sweet[,2])),]

myY.GDD_DTS <- myY[,-2]
myY.GDD_DTS  <- myY.GDD_DTS [-which(is.na(myY.GDD_DTS [,2])),]

#Set the working directory
#setwd(paste(home.dir,"/Pheno+Geno_Data",sep = ""))

#Specify the path where the genotype files are kept
mydatapath <- paste(home.dir,"/Pheno+Geno_Data/", sep="")

#For loop through all of the chromosomess
for(i in 1:10){
  #Set the working directory
  setwd(paste(home.dir,"/Conduct_Kchr/K_chr_Chr_",i,"/",sep=""))
  
  #Read in the kinship matrix and the PCs for the ith chromosome
  myKI <- read.csv("GAPIT.Kin.VanRaden.csv", head = FALSE)
  #head(myKI[,1:7])
  
  myCV <- read.csv("GAPIT.PCA.csv", head = TRUE)
  
  #run GAPIT for sweet corn
  myGAPIT <- GAPIT(
           Y=myY.sweet,  			#This is phenotype data
           #Y=NULL,
           G=NULL,				#This is genotype data,set it to NULL with multiple genotype files
           #GD=myGD,
           #GM=myGM,
           
           file.path=mydatapath,		#The location of genotype files
           #PCA.total=5,
           
           file.from = i,
           file.to = i,
           file.total = 1,
           file.G="AmesUSInbreds_chr", #Common file name (before the numerical part), set it to NULL if for single file
           file.Ext.G="txt.hmp.txt", 
           #file.fragment = 100,
           
          
           
           #SNP.impute= "Major",
           
           #Common file extention (after the numerical part), set it to NULL if for single file
           
           #GDFile=myGDFile,
           #GMFile=myGMFile,
           #GDFileExt=myGDFileExt, 		
           #GMFileExt=myGMFileExt, 		
           
           KI=myKI,				#This is kinship data, set it to NULL in case that geneotype files are used for estimation
           CV=myCV,				#This is the covariate variables of fixed effects, such as population structure
           #Z=myZ,				#This is the customized Z ma
           group.from=2145,		#Was 232	#Lower bound for number of group
           group.to=2145,			#Upper bound for number of group
           group.by=1,				#rang between 1 and number of individuals, smaller the finner optimization 
           #kinship.cluster="average", 			#clustering method: "average", "complete", "ward", "single", "mcquitty", "median","centroid". Example: CA=c("complete","average")
           #kinship.group="Mean",     		#Group kinship tppe:  "Mean", "Max", "Min", "Median". Example: KT=c("Mean","Max")
           #kinship.algorithm="VanRaden", 
           SNP.P3D=TRUE,		#This is the option to use P3D (TRUE) or not (FALSE)
           
           #SNP.impute = "Major",
           SNP.MAF = 0.05,
           cutOff = 0.00,
           #SNP.fraction=0.10
           #Model.selection = TRUE
           
           
           
    )
  
     
  
  
  #run GAPIT for sweet corn
  myGAPIT <- GAPIT(
        Y=myY.GDD_DTS,    		#This is phenotype data
        #Y=NULL,
        G=NULL,				#This is genotype data,set it to NULL with multiple genotype files
        #GD=myGD,
        #GM=myGM,
        
        file.path=mydatapath,		#The location of genotype files
        #PCA.total=5,
        
        file.from = i,
        file.to = i,
        file.total = 1,
        file.G="AmesUSInbreds_chr", #Common file name (before the numerical part), set it to NULL if for single file
        file.Ext.G="txt.hmp.txt", 
        #file.fragment = 100,
        
        
        #SNP.impute= "Major",
        
        #Common file extention (after the numerical part), set it to NULL if for single file
        
        #GDFile=myGDFile,
        #GMFile=myGMFile,
        #GDFileExt=myGDFileExt, 		
        #GMFileExt=myGMFileExt, 		
        
        KI=myKI,				#This is kinship data, set it to NULL in case that geneotype files are used for estimation
        CV=myCV,				#This is the covariate variables of fixed effects, such as population structure
        #Z=myZ,				#This is the customized Z ma
        group.from=1450,		#Was 232	#Lower bound for number of group
        group.to=1450,			#Upper bound for number of group
        #group.by=1,				#rang between 1 and number of individuals, smaller the finner optimization 
        #kinship.cluster="average", 			#clustering method: "average", "complete", "ward", "single", "mcquitty", "median","centroid". Example: CA=c("complete","average")
        #kinship.group="Mean",     		#Group kinship tppe:  "Mean", "Max", "Min", "Median". Example: KT=c("Mean","Max")
        #kinship.algorithm="VanRaden", 
        SNP.P3D=TRUE,		#This is the option to use P3D (TRUE) or not (FALSE)
        
        #SNP.impute = "Major",
        SNP.MAF = 0.05,
        cutOff = 0.00,
        #SNP.fraction=0.10
        #Model.selection = TRUE
    
    
    
  )
  
  #Remove the kinship matrix and PCs from the ith chromosome. This is a safety feature to 
  # make sure that these objects are not inadvertantly used for the (i+1)th chromosome
  rm(myKI)
  rm(myCV)
  
}#End for(i in 1:10)    
     
         
