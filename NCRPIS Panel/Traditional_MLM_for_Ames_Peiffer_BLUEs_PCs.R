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

myY  <- read.table(paste(home.dir, "/Pheno+Geno_Data_Peiffer/Peiffer_et_al_Phenotypes_BLUEs_for_R.txt", sep = ""), head = TRUE)
#head(myY)


#AEL: We do not need to read in the kinship matrix or PCs because we are going to recalculate a separate one for each chromosome. I am also commenting out Lines 64 and 65 so that GAPIT will recalcualte the kinship matrix and PCs from scratch.
myKI <- read.csv(paste(home.dir,"/Pheno+Geno_Data_Peiffer/GAPIT.Kin.VanRaden.csv",sep=""), head = FALSE)
#head(myKI[,1:7])

myCV <- read.csv(paste(home.dir,"/Pheno+Geno_Data_Peiffer/GAPIT.PCA.csv",sep=""), head = TRUE)
#myCV <- myCV[,c(1,ncol(myCV))]
#head(myCV)




#Set the working directory
setwd(paste(home.dir,"/Pheno+Geno_Data_Peiffer",sep = ""))

myY.PHT <- myY[,c(1,2)]
#myY.PHT <- myY.PHT[-which(is.na(myY.PHT[,2])),]
  
  
myY.DTA <- myY[,c(1,3)]
#myY.DTA  <- myY.DTA[-which(is.na(myY.DTA[,2])),]


myY.EHT <- myY[,c(1,4)]
#myY.EHT  <- myY.EHT[-which(is.na(myY.EHT[,2])),]

#NOTE: Things that you have changed:
# file.to and file.from
# group.from = group.to should be 2631


#run GAPIT for plant height
myGAPIT <- GAPIT(
         Y=myY.PHT,  			#This is phenotype data
         #Y=NULL,
         G=NULL,				#This is genotype data,set it to NULL with multiple genotype files
         #GD=myGD,
         #GM=myGM,
         
         #file.path=file.path(mydataPath, subdataPath),		#The location of genotype files
         #PCA.total=5,
         
         file.from = 1,
         file.to = 10,
         file.total = 10,
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
         #group.from=2293,		#Was 232	#Lower bound for number of group
         #group.to=2293,			#Upper bound for number of group
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

   

#run GAPIT for days to anthesis
myGAPIT <- GAPIT(
  Y=myY.DTA,    		#This is phenotype data
  #Y=NULL,
  G=NULL,				#This is genotype data,set it to NULL with multiple genotype files
  #GD=myGD,
  #GM=myGM,
  
  #file.path=file.path(mydataPath, subdataPath),		#The location of genotype files
  #PCA.total=5,
  
  file.from = 1,
  file.to = 10,
  file.total = 10,
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
  #group.from=2145,		#Was 232	#Lower bound for number of group
  #group.to=2145,			#Upper bound for number of group
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

#run GAPIT for ear height
myGAPIT <- GAPIT(
  Y=myY.EHT,      	#This is phenotype data
  #Y=NULL,
  G=NULL,				#This is genotype data,set it to NULL with multiple genotype files
  #GD=myGD,
  #GM=myGM,
  
  #file.path=file.path(mydataPath, subdataPath),		#The location of genotype files
  #PCA.total=5,
  
  file.from = 1,
  file.to = 10,
  file.total = 10,
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
  #group.from=1923,		#Was 232	#Lower bound for number of group
  #group.to=1923,			#Upper bound for number of group
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

 
         