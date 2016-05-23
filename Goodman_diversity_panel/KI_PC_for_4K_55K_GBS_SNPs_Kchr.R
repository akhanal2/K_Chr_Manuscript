#Calculate Kinship and Principal Components for 4K and 55K and GBS SNPs using K_chr 
#######################################################################################
#Step 1: Import files
#######################################################################################
setwd() 
library('MASS')
source("https://bioconductor.org/biocLite.R")
biocLite("multtest")
library(multtest)
library(gplots)
source("http://zzlab.net/GAPIT/emma.txt")
source("GAPIT_Code_from_Internet_20120411_Allelic_Effect.r")
mydata  <- read.delim("SNP55K_maize282_AGPv2_20100513_1_wo_Industry.hmp.txt", head = FALSE) #read in genotype 
myY  <- read.table("mdp_traits.txt", head = TRUE) #read in phenotype (theorectically this step is not necessary 
#because calculations don't depend on phenotype. However, this step is required in order to run GAPIT later)

#######################################################################################
#Step 2: Set result directory and run GAPIT
#######################################################################################
for (i in 1:10){
  #looping through a total of 10 chromosomes 
  
  obs <- which(as.character(mydata[,3])==as.character(i))
  #parse data (extract SNPs from chromosome i)
  
  sub<- mydata[-obs,]
  #create a data file excluding SNPs from chromosome i
  
  subdataPath=paste("Kchr_chrom_", as.character(i), sep="")
  dir.create(subdataPath, showWarnings = FALSE)
  setwd(subdataPath)
  #create a new directory to store kinship and PCs 
                     
  myGAPIT <- GAPIT(
           Y=myY,           #phenotype data
           G=sub,           #genotype data   
           file.path=NULL,		#set the location of genotype files to NULL
           PCA.total = 3,   #take the first 3 PCs
           kinship.algorithm ="Loiselle", #calculate kinship using the Loiselle method
           group.from=278,		#Lower bound for number of group in the phenotype file  
           group.to=278,			#Upper bound for number of group 
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


