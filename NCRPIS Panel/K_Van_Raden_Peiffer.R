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

#myY  <- read.table(paste(home.dir, "/Pheno+Geno_Data/Romay_etal_2013_GenomeBiol_phenotypes-130503_for_R.txt", sep = ""), head = TRUE)
#myY <- myY[,c(1,3)]
#head(myY)


#AEL: We do not need to read in the kinship matrix or PCs because we are going to recalculate a separate one for each chromosome. I am also commenting out Lines 64 and 65 so that GAPIT will recalcualte the kinship matrix and PCs from scratch.
#myKI <- read.csv(paste(mydataPath,"Comprehensive_K.csv",sep=""), head = FALSE)
#head(myKI[,1:7])

#myCV <- read.csv(paste(mydataPath,"GAPIT.PCA.csv",sep=""), head = TRUE)
#myCV <- myCV[,c(1,ncol(myCV))]
#head(myCV)


#Set the working directory
setwd(paste(home.dir,"/Pheno+Geno_Data_Peiffer", sep = ""))

#run GAPIT (calculate K using Loiselle)
myGAPIT <- GAPIT(
         #Y=myY,  			#This is phenotype data
         Y=NULL,
         G=NULL,				#This is genotype data,set it to NULL with multiple genotype files
         #GD=myGD,
         #GM=myGM,
         
         #file.path=file.path(mydataPath, subdataPath),		#The location of genotype files
         PCA.total=5,
         
         file.from = 1,
         file.to = 10,
         file.total = 10,
         file.G="AmesUSInbreds_chr", #Common file name (before the numerical part), set it to NULL if for single file
         file.Ext.G="txt.hmp.txt", 
         
        
         
         #SNP.impute= "Major",
         
         #Common file extention (after the numerical part), set it to NULL if for single file
         
         #GDFile=myGDFile,
         #GMFile=myGMFile,
         #GDFileExt=myGDFileExt, 		
         #GMFileExt=myGMFileExt, 		
         
         #KI=myKI,				#This is kinship data, set it to NULL in case that geneotype files are used for estimation
         #CV=myCV,				#This is the covariate variables of fixed effects, such as population structure
         #Z=myZ,				#This is the customized Z ma
         group.from=10,		#Was 232	#Lower bound for number of group
         group.to=10,			#Upper bound for number of group
         group.by=1,				#rang between 1 and number of individuals, smaller the finner optimization 
         kinship.cluster="average", 			#clustering method: "average", "complete", "ward", "single", "mcquitty", "median","centroid". Example: CA=c("complete","average")
         kinship.group="Mean",     		#Group kinship tppe:  "Mean", "Max", "Min", "Median". Example: KT=c("Mean","Max")
         kinship.algorithm="VanRaden", 
         SNP.P3D=TRUE,		#This is the option to use P3D (TRUE) or not (FALSE)
         
         SNP.impute = "Major",
         SNP.MAF = 0.05,
         cutOff = 0.00,
         SNP.fraction=0.10
         #Model.selection = TRUE
         
         
         
       )

  

########################Calcualte a separate kinship matrix within each chromosome

for(i in 1:10){
  setwd(paste(home.dir,"/Conduct_K_chr_Peiffer/K_chr_Chr_",i,"/", sep = ""))
  
  #run GAPIT (calculate K using Loiselle)
  myGAPIT <- GAPIT(
    #Y=myY,    		#This is phenotype data
    Y=NULL,
    G=NULL,				#This is genotype data,set it to NULL with multiple genotype files
    #GD=myGD,
    #GM=myGM,
    
    #file.path=file.path(mydataPath, subdataPath),		#The location of genotype files
    PCA.total=5,
    
    file.from = 1,
    file.to = 9,
    file.total = 9,
    file.G="AmesUSInbreds_Kchr", #Common file name (before the numerical part), set it to NULL if for single file
    file.Ext.G="txt.hmp.txt", 
    
    
    
    #SNP.impute= "Major",
    
    #Common file extention (after the numerical part), set it to NULL if for single file
    
    #GDFile=myGDFile,
    #GMFile=myGMFile,
    #GDFileExt=myGDFileExt, 		
    #GMFileExt=myGMFileExt, 		
    
    #KI=myKI,				#This is kinship data, set it to NULL in case that geneotype files are used for estimation
    #CV=myCV,				#This is the covariate variables of fixed effects, such as population structure
    #Z=myZ,				#This is the customized Z ma
    group.from=10,		#Was 232	#Lower bound for number of group
    group.to=10,			#Upper bound for number of group
    group.by=1,				#rang between 1 and number of individuals, smaller the finner optimization 
    kinship.cluster="average", 			#clustering method: "average", "complete", "ward", "single", "mcquitty", "median","centroid". Example: CA=c("complete","average")
    kinship.group="Mean",     		#Group kinship tppe:  "Mean", "Max", "Min", "Median". Example: KT=c("Mean","Max")
    kinship.algorithm="VanRaden", 
    SNP.P3D=TRUE,		#This is the option to use P3D (TRUE) or not (FALSE)
    
    SNP.impute = "Major",
    SNP.MAF = 0.05,
    cutOff = 0.00,
    SNP.fraction=0.10
    #Model.selection = TRUE
    
    
    
  )

}#end  for(i in 1:10) 
     
         