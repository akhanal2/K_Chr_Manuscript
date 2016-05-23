#Step 1: Set data directory and import files
#######################################################################################
setwd("/Users/adminuser/Desktop/Work/Graduate_Students/Angela_(Hsaio-Han)_Chen/Data_for_Angela/Ames_Data")
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

#AEL: We do not need to read in the kinship matrix or PCs because we are going to recalculate a separate one for each chromosome. I am also commenting out Lines 64 and 65 so that GAPIT will recalcualte the kinship matrix and PCs from scratch.
myKI <- read.csv(paste(home.dir,"/Pheno+Geno_Data/GAPIT.Kin.VanRaden.csv",sep=""), head = FALSE)
#head(myKI[,1:7])

myCV <- read.csv(paste(home.dir,"/Pheno+Geno_Data/GAPIT.PCA.csv",sep=""), head = TRUE)
#myCV <- myCV[,c(1,ncol(myCV))]
#head(myCV)




#Set the working directory
setwd(paste(home.dir,"/Pheno+Geno_Data",sep = ""))

#myY.sweet <- myY[,-3]
#myY.GDD_DTS <- myY[,-2]

#5 PCs
myGAPIT <- GAPIT(
  Y=myY.sweet,        #This is phenotype data
  #Y=NULL,
  G=NULL,    		#This is genotype data,set it to NULL with multiple genotype files
  
  file.from = 1,
  file.to = 1,
  file.total = 1,
  file.G="AmesUSInbreds_chr_tmp", #Common file name (before the numerical part), set it to NULL if for single file
  file.Ext.G="txt.hmp.txt", 
  
  
  KI=myKI,				#This is kinship data, set it to NULL in case that geneotype files are used for estimation
  CV=myCV,				#This is the covariate variables of fixed effects, such as population structure
  #Z=myZ,				#This is the customized Z ma
  group.from=2145,		#Was 232	#Lower bound for number of group
  group.to=2145,			#Upper bound for number of group
  
  SNP.P3D=TRUE,	
  
  SNP.MAF = 0.05,
  cutOff = 0.00, 
)


         