############Summary Table for Results#############
#setwd("C:/Users/chen398/Desktop/New_Results")

setwd()
##**********Number of Significant Association at 5%/10% FDR**********

#$K_chr Approach$

  #-----Flowering Time-------
dataPath="Conduct_Kchr_Romay/"
traitnames <- c("GDD_DTS","Sweet")
total.count.1<-0
total.count.2<-0
for (i in traitnames){
  results <-read.csv(paste(dataPath,"GAPIT.",i,".GWAS.Results.csv",sep=""),head=TRUE)
  sig.pvalue.count.1 <- nrow(results[which(results[,11]<=0.05),]) ##counts of significant SNPs at 5% FDR
  sig.pvalue.count.2 <- nrow(results[which(results[,11]<=0.10),]) ##counts of significant SNPs at 10% FDR
  total.count.1 <- sig.pvalue.count.1+total.count.1
  total.count.2 <- sig.pvalue.count.2+total.count.2
  write.table(c(paste("Number of significant Associations at 5% FDR: ", total.count.1, sep = ""),
                paste("Number of significant Associations at 10% FDR: ", total.count.2, sep = "")), 
                paste("Number_of_Sig_Associations.K_Chr.", i,".txt", sep = ""), quote = FALSE, 
              sep = "\t", row.names = FALSE,col.names = TRUE)
  }


#sig.pvalue.SNP <- results$SNP[which(results$P.value<= 0.05)] ##names of significant SNPs at 5% FDR

#$Traditional Q+K$

  #-----Flowering Time-------
dataPath="Pheno+Geno_Data_Romay/"
traitnames <- c("GDD_DTS","Sweet")
total.count.1<-0
total.count.2<-0
for (i in traitnames){
  results <-read.csv(paste(dataPath,"GAPIT.",i,".GWAS.Results.csv",sep=""),head=TRUE)
  sig.pvalue.count.1 <- nrow(results[which(results[,10]<=0.05),]) ##counts of significant SNPs at 5% FDR
  sig.pvalue.count.2 <- nrow(results[which(results[,10]<=0.10),]) ##counts of significant SNPs at 10% FDR
  total.count.1 <- sig.pvalue.count.1+total.count.1
  total.count.2 <- sig.pvalue.count.2+total.count.2
  write.table(c(paste("Number of significant Associations at 5% FDR: ", total.count.1, sep = ""),
                paste("Number of significant Associations at 10% FDR: ", total.count.2, sep = "")), 
              paste("Number_of_Sig_Associations.Trad.", i,".txt", sep = ""), quote = FALSE, 
              sep = "\t", row.names = FALSE,col.names = TRUE)
}

#####################################################################################################
##**********Number of Significant Association at 10% FDR in Novel Genomic Regions for K_chr**********


#-----Flowering Time-------
dataPath.1="Conduct_Kchr_Romay/"
dataPath.2="Pheno+Geno_Data_Romay/"
traitnames <-c("GDD_DTS","Sweet")
for (i in traitnames){
  
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.",i,".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.",i,".GWAS.Results.csv",sep=""),head=TRUE)
  sig.pvalue.Kchr <-results.Kchr[which(results.Kchr[,11]<=0.10),]
  sig.pvalue.class <-results.class[which(results.class[,10]<=0.10),]
  novel.genom.chrom <- sig.pvalue.Kchr$Chromosome
  novel.genom.upstr <- sig.pvalue.Kchr$Position+250000
  novel.genom.dwstr <- sig.pvalue.Kchr$Position-250000

  count <-0
  
  for(j in 1:nrow(sig.pvalue.Kchr)){  
    if(nrow(sig.pvalue.Kchr)==0)break
    sig.pvalue.class.in.novel.genom.region<- sig.pvalue.class[which((sig.pvalue.class$Chromosome==novel.genom.chrom[j])
                                                                    &(sig.pvalue.class$Position>=novel.genom.dwstr[j])
                                                                    &(sig.pvalue.class$Position<=novel.genom.upstr[j])),]
    #AEL: I changed this a little bit. Please let me know if you disagree.
    if (nrow(sig.pvalue.class.in.novel.genom.region)==0){ #I chnaged this form "==" to "!=" because my understanding is that
      # novel genomic regions will be flagged in "true.novel.genom.region" and "cum.novel.genom.region". When novel 
      # genomic regions are identified by K_chr, the number of rows of "sig.pvalue.class.in.novel.genom.region"
      # will be greater than 0.
      if(count == 0){
        true.novel.genom.region <- sig.pvalue.Kchr[j,]
      }else{
        true.novel.genom.region <- rbind(true.novel.genom.region, sig.pvalue.Kchr[j,])
      }
      count <-count+1
      cum.novel.genom.region <- nrow(true.novel.genom.region)
    }
    
    
  } #end for(j in 1:nrow(sig.pvalue.Kchr))
  
  if(count!= 0){
    
    #Export "true.novel.genomic.region" to a tab delimited text file
    write.table(true.novel.genom.region, paste("Novel.genomic.regions.detected.for.FT.", i,".txt", sep = ""), quote = FALSE, 
                sep = "\t", row.names = FALSE,col.names = TRUE)
    
    #Export "cum.novel.genom.region" to a tab delmited text file
    write.table(cum.novel.genom.region, paste("Number.of.Novel.genomic.regions.detected.for.FT.", i,".txt", sep = ""), quote = FALSE, 
                sep = "\t", row.names = FALSE,col.names = TRUE)
  }
  
}#end for (i in traitnames)




##**********Number of Significant Association at 10% FDR in Novel Genomic Regions for the traditional MLMr**********



#-----Flowering Time-------
dataPath.1="Conduct_Kchr_Romay/"
dataPath.2="Pheno+Geno_Data_Romay/"
traitnames <- c("GDD_DTS","Sweet")
for (i in traitnames){
  
  results.Kchr <-read.csv(paste(dataPath.1,"GAPIT.",i,".GWAS.Results.csv",sep=""),head=TRUE)
  results.class <-read.csv(paste(dataPath.2,"GAPIT.",i,".GWAS.Results.csv",sep=""),head=TRUE)
  sig.pvalue.Kchr <-results.Kchr[which(results.Kchr[,11]<=0.10),]
  sig.pvalue.class <-results.class[which(results.class[,10]<=0.10),]
  novel.genom.chrom <- sig.pvalue.Kchr$Chromosome
  novel.genom.upstr <- sig.pvalue.Kchr$Position+250000
  novel.genom.dwstr <- sig.pvalue.Kchr$Position-250000
  
  count <-0
  
  for(j in 1:nrow(sig.pvalue.class)){  
    if(nrow(sig.pvalue.class)==0)break
    sig.pvalue.Kchr.in.novel.genom.region<- sig.pvalue.Kchr[which((sig.pvalue.Kchr$Chromosome==novel.genom.chrom[j])
                                                                  &(sig.pvalue.Kchr$Position>=novel.genom.dwstr[j])
                                                                  &(sig.pvalue.Kchr$Position<=novel.genom.upstr[j])),]
    #AEL: I changed this a little bit. Please let me know if you disagree.
    if (nrow(sig.pvalue.Kchr.in.novel.genom.region)
        ==0){ #I chnaged this form "==" to "!=" because my understanding is that
      # novel genomic regions will be flagged in "true.novel.genom.region" and "cum.novel.genom.region". When novel 
      # genomic regions are identified by K_chr, the number of rows of "sig.pvalue.class.in.novel.genom.region"
      # will be greater than 0.
      if(count == 0){
        true.novel.genom.region <- sig.pvalue.class[j,]
      }else{
        true.novel.genom.region <- rbind(true.novel.genom.region, sig.pvalue.class[j,])
      }
      count <-count+1
      cum.novel.genom.region <- nrow(true.novel.genom.region)
    }
    
    
  } #end for(j in 1:nrow(sig.pvalue.Kchr))
  
  
  if(count!= 0){
    
    #Export "true.novel.genomic.region" to a tab delimited text file
    write.table(true.novel.genom.region, paste("Novel.genomic.regions.Trad.detected.", i,".txt", sep = ""), quote = FALSE, 
                sep = "\t", row.names = FALSE,col.names = TRUE)
    
    #Export "cum.novel.genom.region" to a tab delmited text file
    write.table(cum.novel.genom.region, paste("Number.of.Novel.genomic.regions.Trad.detected.", i,".txt", sep = ""), quote = FALSE, 
                sep = "\t", row.names = FALSE,col.names = TRUE)
  }
  
}#end for (i in traitnames)




