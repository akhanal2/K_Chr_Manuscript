**KI_PC_for_4K_55K_GBS_SNPs_Kchr.r** produces a kinship matrix and 3 principal components using the K_Chr method. It reads in 
*"GAPIT_Code_from_Internet_20120411_Allelic_Effect.r"* and 

**_Genotype data_** from *"SNP55K_maize282_AGPv2_20100513_1_wo_Industry.hmp.txt"*

**_Phenotype data_** from one of the following: 

Tocochromanols: *"BLUPs_NO_Outliers_282_Transformed_all.txt"* (number of group: 252)

Carotenoids: *"BLUPs_Transformed_all35.csv"* (number of group: 202)

Flowering time: *"mdp_traits.csv"* (number of group: 278)
> Note that whichever phenotype file you choose doesn't matter as the calculation of kinship and PCs don't depend on phenotype information. However, this step can't be omitted, otherwise code won't run.  


**KI_PC_for_4K_55K_GBS_SNPs_tradMLM.r** produces a kinship matrix and 3 principal components using the traditional MLM method. It reads in the same information described previously for *"KI_PC_for_4K_55K_GBS_SNPs_Kchr.r"*


**GWAS_from_4K_55K_GBS_SNPs.r** takes in the kinships and PCs generated previously and produce GWAS results from different **_genotyping sources_**:

4K SNPs: *"4K_SNPsmdp_genotype_test1.hmp.txt"*

55K SNPs: *"SNP55K_maize282_AGPv2_20100513_1hmp.txt"*

GBS SNPs: 10 files in total *"282_20111217_scv10mF8maf002_mgs_E1pLD5.imp95_1024.c1.hmp.rar"*......*"282_20111217_scv10mF8maf002_mgs_E1pLD5.imp95_1024.c10.hmp.rar"* 

and the aforementioned **_phenotype data_**: Tocochromanols, Carotenoids, and Flowering time.   

> Note that GBS files contain taxa names formatted differently from 4K and 55K. To make life easier, please use the following **_phenotype_** files formatted specifically for GBS:

> Tocochromanols: *"BLUPs_NO_Outliers_282_Transformed_all_GBS_Names.txt"*

> Carotenoids: *"BLUPs_Transformed_all35_GBSNames.csv"*

> Flowering time: *"mdp_traits_GBS.csv"* 

This code requires *"mt.rawp2adjp.r"* and *"GAPIT_Code_from_Internet_20120411_Allelic_Effect.r"* to run. You can generate traditional MLM GWAS results with this code using kinships and PCs generated from *"KI_PC_for_4K_55K_GBS_SNPs_tradMLM.r"*.


**Combine_GWAS_from_4K_55K_GBS.r** combines all the GWAS results generated previously. This file reads in   *"GAPIT_with_BIC_20120210.r"* 


**K_Chr_Performance.r** finds novel genomic regions defined in the paper 


**Wilcoxon.r** compares significant marker (from specific genomic regions) distributions of K_chr and traditional MLM 


**Box Plots.r** creates boxplots for significant markers identified using K_chr vs traditional MLM


**Create_Manhattan_Plot.R** creates manhattan plots for tocochromanols and carotenoids. Supporting files are in *"Tocochromanol_Traits_in_Novel_Genomic_Regions.rar"* and *"Carotenoid_Traits_in_Novel_Genomic_Regions.rar"*. 

