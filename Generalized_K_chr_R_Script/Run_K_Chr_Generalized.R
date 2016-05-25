home.dir <- getwd()


library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

#Read in "K_Chr_Generalized.R". This file, which contains the "GAPIT.Kchr" funciton (which runs the K_chr model),
# should be in the same directory as the directory indicated in line 1 of this script.
source("K_Chr_Generalized.R")


#Create a subdirectory where all input files are stored, and read in the phenotypic data
mydatapath <- paste(home.dir,"/GAPIT_Demo_Data/", sep="")
myY  <- read.table(paste(mydatapath, "/mdp_traits.txt", sep = ""), head = TRUE)


#Run the K_chr model using the GAPIT.Kchr() funciton
myGAPIT.Kchr <- GAPIT.Kchr(mydatapath = mydatapath, 
                           myY  = myY ,
                           myfile.G = "mdp_genotype_chr", #Common file name (before the numerical part), set it to NULL if for single file
                           myfile.Ext.G = "hmp.txt",
                           myfile.from = 1,
                           myfile.to = 10,
                           mykinship.algorithm = "Loiselle",
                           mygroup.from = 281,
                           mygroup.to = 281,
                           myPCA.total = 5,
                           mySNP.impute = "Major",
                           mySNP.MAF = 0.05,
                           mycutOff = 0.00,
                           mySNP.fraction = 1,
                           myModel.selection = TRUE,
                           mySNP.P3D = TRUE,
                           home.dir = home.dir
)






     
         
