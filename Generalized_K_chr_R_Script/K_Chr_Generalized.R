
#########################Function created by Alex Lipka and Angela Chen, 4/6/2016
####Required input files:
##### Genotypic data in hapmap format; with a separate file for each chromosome
##### Phenotypic data
##### Please see the GAPIT user manual for a description of how to format the input file
GAPIT.Kchr <- function(mydatapath = NULL, home.dir = NULL, myY  = NULL, myfile.G = NULL, myfile.Ext.G = NULL, myfile.from = NULL, myfile.to = NULL,
                       mykinship.algorithm = "VanRaden", mygroup.from = 0, mygroup.to = 10000, myPCA.total = 0, mySNP.impute = "Middle",
                       mySNP.MAF = 0.00, mycutOff = 0.01, mySNP.fraction = 1, myModel.selection = FALSE, mySNP.P3D = TRUE){
                       
#############    
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
    if(myPCA.total != 0){
      myCV <- read.csv("GAPIT.PCA.csv", head = TRUE)
      myCV <- myCV[,c(1:(myPCA.total+1))]
    }
    else{
      myCV = NULL
    }

  
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
}#end GAPIT.Kchr
###########################
