setwd("/Users/adminuser/Desktop/Work/Graduate_Students/Angela_(Hsaio-Han)_Chen/Genearlized_K_chr_Script/")
home.dir <- getwd()


library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

#Create a subdirectory where all input files are stored, and read in the phenotypic data
mydatapath <- paste(home.dir,"/GAPIT_Demo_Data/", sep="")
myY  <- read.table(paste(mydatapath, "/mdp_traits.txt", sep = ""), head = TRUE)

#########################Function created by Alex Lipka, 4/6/2016
####Required input files:
##### Genotypic data in hapmap format; with a separate file for each chromosome
##### Phenotypic data
##### Please see the GAPIT user manual for a description of how to format the input file
mydatapath <- mydatapath 
myY  <- myY 
myfile.G <- "mdp_genotype_chr" #Common file name (before the numerical part), set it to NULL if for single file
myfile.Ext.G <- "hmp.txt"
myfile.from <- 1
myfile.to <- 10
mykinship.algorithm <- "Loiselle"
mygroup.from <- 281
mygroup.to <- 281
mygroup.by <- 1
myPCA.total <- 5
mySNP.impute <- "Major"
mySNP.MAF <- 0.05
mycutOff <- 0.00
mySNP.fraction <- 1
myModel.selection <- TRUE
mySNP.P3D <- TRUE

#### Eventually include all possible input parameters that the GAPIT function has


for(i in myfile.from:myfile.to){ #For loop through each chromosome
  #Create a directory
  dir.create(paste("K_chr_Chr", i, "", sep = ""))
  new.path.dirname <- paste(home.dir, "/K_chr_Chr", i,sep = "")
  
  #Put the all of the genotype files except for the one being tested in this directory
  count <- 1
  sequence.of.numbers.for.below.loop <- seq(myfile.from:myfile.to)[-i]
  for(j in sequence.of.numbers.for.below.loop){
    file.copy(paste(mydatapath, myfile.G,j,".",myfile.Ext.G, sep=""),   paste(new.path.dirname, "/", myfile.G,count,".",myfile.Ext.G, sep=""))
    count <- count+1
  }

  
  #Calculate the appropriate kinship matrix and PCs
  setwd(new.path.dirname)
  
  #run GAPIT (calculate K using Loiselle)
  myGAPIT <- GAPIT(
    Y=NULL,

    PCA.total=myPCA.total,
    
    file.from = 1,
    file.to = (count-1),
    file.total = length(sequence.of.numbers.for.below.loop),
    file.G=myfile.G, #Common file name (before the numerical part), set it to NULL if for single file
    file.Ext.G=myfile.Ext.G, 
    
 
    group.from=10,		#Was 232	#Lower bound for number of group
    group.to=10,			#Upper bound for number of group
    group.by=1,				#range between 1 and number of individuals, smaller the finner optimization 
    kinship.cluster="average", 			#clustering method: "average", "complete", "ward", "single", "mcquitty", "median","centroid". Example: CA=c("complete","average")
    kinship.group="Mean",     		#Group kinship tppe:  "Mean", "Max", "Min", "Median". Example: KT=c("Mean","Max")
    kinship.algorithm=mykinship.algorithm, 
    SNP.P3D=TRUE,		#This is the option to use P3D (TRUE) or not (FALSE)
    
    SNP.impute = mySNP.impute,
    SNP.MAF = mySNP.MAF,
    cutOff = mycutOff,
    SNP.fraction=mySNP.fraction,  
    Geno.View.output=FALSE
  )
  
  print(paste("------------------------------- Chromosome ", i, ": Just finished calcualting separate kinship matrices, now performing association analysis on this chromosome... -----------------------------  ", sep = ""))
  
  #Read in the appropriate kinship matrix and PCs that were just calcualted
  myKI <- read.csv(paste("GAPIT.Kin.",mykinship.algorithm,".csv", sep = ""), head = FALSE)
  myCV <- read.csv("GAPIT.PCA.csv", head = TRUE)
  myCV <- myCV[,c(1:(myPCA.total+1))]

  #Run GWAS on that chromosome
  
  myGAPIT <- GAPIT(
    Y=myY,    		#This is phenotype data
    G=NULL,				#This is genotype data,set it to NULL with multiple genotype files   
    file.path=mydatapath,		#The location of genotype files
    
    file.from = i,
    file.to = i,
    file.total = 1,
    file.G=myfile.G, #Common file name (before the numerical part), set it to NULL if for single file
    file.Ext.G=myfile.Ext.G, 
   
  
  
    KI=myKI,				#This is kinship data, set it to NULL in case that geneotype files are used for estimation
    CV=myCV,				#This is the covariate variables of fixed effects, such as population structure
    #Z=myZ,				#This is the customized Z ma
    group.from=mygroup.from ,		#Was 232	#Lower bound for number of group
    group.to=mygroup.from,			#Upper bound for number of group
    #group.by=1,
  
    SNP.P3D=mySNP.P3D,		#This is the option to use P3D (TRUE) or not (FALSE)
    
    SNP.impute = mySNP.impute,
    SNP.MAF = mySNP.MAF,
    cutOff = mycutOff,
    Geno.View.output=FALSE,

    Model.selection = myModel.selection
    
  )
  
  setwd(home.dir)

}#End for(i in myfile.from:myfile.to)


#Run that function that redoes the Benjamini-Hochberg (B-H) FDR-controlling procedure for the genome-wide SNPs, and remakes
# the genome-wide Manhattan plots, etc.

name.of.traits <- colnames(myY)[-1]
#Read in the separate results from each chromosome
for(i in name.of.traits){
  
  #This count object will be used in the for loop to help efficiently append GWAS results
  # into one object
  count <- 0
  for(j in myfile.from:myfile.to){#For loop through the chromosomes
    
    #Read in the GWAS results
    GWAS.results.this.chr  <- read.csv(paste(home.dir,"/K_chr_Chr",j,"/GAPIT..",i,".GWAS.Results.csv",sep=""), head = TRUE)
    
    
    #Append them to the previous set of GWAS results
    if(count == 0){
      GWAS.results <- GWAS.results.this.chr
    }else{
      GWAS.results <- rbind(GWAS.results, GWAS.results.this.chr)
    }#end if(count == 0)
    
    count <- count+1
    
  } #End for(j in 1:10)
  
  
  #Run the B-H procedure on the combined data, and append the FDR-adjusted P-values to the GWAS.Results
  Conduct.FDR <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = GWAS.results, 
                                                                    FDR.Rate = 0.05, FDR.Procedure = "BH")
  
  GWAS.Results.with.FDR <- Conduct.FDR$PWIP
  
  #Make new Manhattan plots with the combined results
  GAPIT.Manhattan(GI.MP = GWAS.Results.with.FDR[,2:4], name.of.trait = paste(i,".K_chr",sep = ""), 
                  DPP=50000, plot.type = "Genomewise",cutOff=0.00)
  
  GAPIT.Manhattan(GI.MP = GWAS.Results.with.FDR[,2:4], name.of.trait = paste(i,".K_chr",sep = ""), 
                  DPP=50000, plot.type = "Chromosomewise",cutOff=0.00)
  
  colnames(GWAS.Results.with.FDR)[(ncol(GWAS.Results.with.FDR)-1)] = "Chromosome.wide.FDR.adj.P"
  
  colnames(GWAS.Results.with.FDR)[ncol(GWAS.Results.with.FDR)] = "Genomewide.wide.FDR.adj.P"
  
  write.table(GWAS.Results.with.FDR, paste("GAPIT.", i,".K_chr.GWAS.Results.csv", sep = ""), 
              quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
  
  
  rm(GWAS.results.this.chr)
  rm(GWAS.results)
  rm(GWAS.Results.with.FDR)
  rm(Conduct.FDR)
} #End for(i in traits.total)




















#Old code, which can be discarded
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

myY  <- read.table(paste(home.dir, "/Pheno+Geno_Data_Peiffer/Peiffer_et_al_Phenotypes_BLUPs_for_R.txt", sep = ""), head = TRUE)
#head(myY)

#Set the working directory
setwd(paste(home.dir,"/Pheno+Geno_Data_Peiffer",sep = ""))

myY.PHT <- myY[,c(1,2)]
#myY.PHT <- myY.PHT[-which(is.na(myY.PHT[,2])),]


myY.DTA <- myY[,c(1,3)]
#myY.DTA  <- myY.DTA[-which(is.na(myY.DTA[,2])),]


myY.EHT <- myY[,c(1,4)]
#myY.EHT  <- myY.EHT[-which(is.na(myY.EHT[,2])),]

#Set the working directory
#setwd(paste(home.dir,"/Pheno+Geno_Data",sep = ""))

#Specify the path where the genotype files are kept
mydatapath <- paste(home.dir,"/Pheno+Geno_Data_Peiffer/", sep="")

#For loop through all of the chromosomess
for(i in 1:10){
  #Set the working directory
  setwd(paste(home.dir,"/Conduct_K_chr_Peiffer/K_chr_Chr_",i,"/",sep=""))
  
  #Read in the kinship matrix and the PCs for the ith chromosome
  myKI <- read.csv("GAPIT.Kin.VanRaden.csv", head = FALSE)
  #head(myKI[,1:7])
  
  #Conduct spectral decomposition on the kinship matrix
  spectral.decomp.of.myKI = eigen(myKI[,-1], symmetric=TRUE, only.values = FALSE)
  first.six.eigenvectors.of.myKI <- spectral.decomp.of.myKI$vectors[,c(1:6)]
  
  
  
  myCV <- data.frame(as.character(myKI[,1]), first.six.eigenvectors.of.myKI)
  
  #run GAPIT for plant height
  myGAPIT <- GAPIT(
           Y=myY.PHT,  			#This is phenotype data
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
           group.from=1903,		#Was 232	#Lower bound for number of group
           group.to=1903,			#Upper bound for number of group
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
  
  
  #run GAPIT for DTA
  myGAPIT <- GAPIT(
    Y=myY.DTA,    		#This is phenotype data
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
    group.from=1423,		#Was 232	#Lower bound for number of group
    group.to=1423,			#Upper bound for number of group
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
  

  #run GAPIT for Ear Height
  myGAPIT <- GAPIT(
    Y=myY.EHT,      	#This is phenotype data
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
    group.from=1923,		#Was 232	#Lower bound for number of group
    group.to=1923,			#Upper bound for number of group
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
  
  
  
  
  #Remove the kinship matrix and PCs from the ith chromosome. This is a safety feature to 
  # make sure that these objects are not inadvertantly used for the (i+1)th chromosome
  rm(myKI)
  rm(myCV)
  
}#End for(i in 1:10)    
     
         