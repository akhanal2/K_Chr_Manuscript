# GAPIT - Genome Association and Prediction Integrated Tool
# Designed by Zhiwu Zhang
# Written by Zhiwu Zhang, Alex Lipka and Feng Tian 
# Last update: August 10, 2011 

#Step 0: Set directory to load GAPIT source files (this step is omited for using package)
#######################################################################################
library('MASS')
library(multtest)
library(gplots)

source("http://zzlab.net/GAPIT/emma.txt")
source("GAPIT_with_BIC_20120210.r")
#source("http://www.maizegenetics.net/images/stories/bioinformatics/GAPIT/gapit_functions.txt")
traits.total <- c("ALL_PHT_BLUE","ALL_DTA_BLUE", "ALL_EHT_BLUE")
number.of.chromosomes <- 10

#Step 1: Set data directory and import files
#######################################################################################
#mydataPath.GBS="D:\\Paper_with_Brenda\\GWAS_Tocochromanols\\3PC+K_282_BLUPs_Only_GBS\\"
#mydataPath.55K="D:\\Paper_with_Brenda\\GWAS_Tocochromanols\\3PC+K_282_BLUPs_Only_55K\\"
#mydataPath.4K ="D:\\Paper_with_Brenda\\GWAS_Tocochromanols\\3PC+K_282_BLUPs_Only_4K\\"

#Transformed <- TRUE

#if(Transformed) trait <- "traittrans"
#if(!Transformed) trait <- "trait"

#Step 2: Set the result directory to where the combined data will go
#######################################################################################
setwd("/Users/adminuser/Desktop/Work/Graduate_Students/Angela_(Hsaio-Han)_Chen/Data_for_Angela/Ames_Data/Conduct_K_chr_Peiffer")
home.dir <- getwd()

for(i in traits.total){
  
  #This count object will be used in the for loop to help efficiently append GWAS results
  # into one object
  count <- 0
  for(j in 1:number.of.chromosomes){#For loop through the chromosomes
  
    #Read in the GWAS results
    GWAS.results.this.chr  <- read.csv(paste(home.dir,"/K_chr_Chr_",j,"/GAPIT.",i,".GWAS.Results.csv",sep=""), head = TRUE)
  
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
  GAPIT.Manhattan(GI.MP = GWAS.Results.with.FDR[,2:4], name.of.trait = i, 
                  DPP=50000, plot.type = "Genomewise",cutOff=0.00)
  
  GAPIT.Manhattan(GI.MP = GWAS.Results.with.FDR[,2:4], name.of.trait = i, 
                  DPP=50000, plot.type = "Chromosomewise",cutOff=0.00)
  
  
  write.table(GWAS.Results.with.FDR, paste("GAPIT.", i,".GWAS.Results.csv", sep = ""), 
              quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
  

  rm(GWAS.results.this.chr)
  rm(GWAS.results)
  rm(GWAS.Results.with.FDR)
  rm(Conduct.FDR)
} #End for(i in traits.total)




#
