#Compile GWAS results from 4K and 55K and GBS SNPs
#######################################################################################
#Step 1: Import files
#######################################################################################
setwd()
library('MASS')
library(multtest)
library(gplots)
source("http://zzlab.net/GAPIT/emma.txt")
source("GAPIT_with_BIC_20120210.r")
traits.total= 
  #traits.total=c(1:20) for tocochromanol traits
  #traits.total=c(1:4,6,7,10,13,14,16,17,21,25,28,30) for carotenoid traits  
  #traits.total=c(dpoll,EarHT,EarDia) for flowering time traits
#######################################################################################
#Step 2: Set result directory 
#######################################################################################
subdataPath="Combined"
dir.create(subdataPath, showWarnings = FALSE)
setwd(subdataPath)
#create a directory to store results


for(i in traits.total){  
  count <- 0 #initiate count value to 0
  for(j in 1:10){  #loop through 10 chromosomes to create a dataset with all chromosomes (1-10) results 
  
    GWAS.4K  <- read.csv(paste("Kchr_GWAS_4K_chrom_",j,"/GAPIT.",trait,i ,".GWAS.Results.csv",sep=""), head = TRUE)
    #read in 4K SNPs GWAS results from chrom j 
    
    GWAS.55K  <- read.csv(paste("Kchr_GWAS_55K_chrom_",j,"/GAPIT.",trait,i ,".GWAS.Results.csv",sep=""), head = TRUE)
    #read in 55K SNPs GWAS results from chrom j 
    
    GWAS.GBS  <- read.csv(paste(mydataPath,"Kchr_GWAS_GBS_chrom_",j,"/GAPIT.",trait,i ,".GWAS.Results.csv",sep=""), head = TRUE)
    #read in GBS SNPs GWAS results from chrom j 
    
    results.this.chr <- rbind(GWAS.GBS, GWAS.55K, GWAS.4K) 
    #combine all the GWAS results from chrom j
    
    if(count == 0){
      GWAS.Results <- results.this.chr
      }else{
      GWAS.Results <- rbind(GWAS.Results, results.this.chr)
      } #create a dataset with all chromosomes (1-10) results 
    count <- count+1
  } #end loop for for(j in 1:10)
  
  
  GWAS.Results <- GWAS.Results[,-ncol(GWAS.Results)]
  #Remove the last column, which is the FDR adjusted P-values
  
  Conduct.FDR <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = GWAS.Results, 
                           FDR.Rate = 0.05, FDR.Procedure = "BH")
  GWAS.Results.with.FDR <- Conduct.FDR$PWIP
  #Run the B-H procedure on the combined data, and append the FDR-adjusted P-values to the GWAS.Results
  
  GAPIT.Manhattan(GI.MP = GWAS.Results.with.FDR[,2:4], name.of.trait = paste(trait, i, sep = ""), 
                  DPP=50000, plot.type = "Genomewise",cutOff=0.00)
  GAPIT.Manhattan(GI.MP = GWAS.Results.with.FDR[,2:4], name.of.trait = paste(trait, i, sep = ""), 
                  DPP=50000, plot.type = "Chromosomewise",cutOff=0.00)
  #Make new Manhattan plots with the combined results
  
  write.table(GWAS.Results.with.FDR, paste("GAPIT.", trait, i, ".GWAS.Results.csv", sep = ""), 
                 quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
  #create a new csv file with the combined results
  
  rm(GWAS.GBS)
  rm(GWAS.55K)
  rm(GWAS.4K)
  rm(results.this.chr)
  rm(GWAS.Results)
  rm(GWAS.Results.with.FDR)
  rm(Conduct.FDR)

} #End loop for for(i in traits.total)




