# K_Chr_Manuscript
 "The Use of Targeted Marker Subsets to Account for Population Structure and Relatedness in Genome-Wide Association Studies of Maize (Zea mays L.)" manuscript by Angela H. Chen and Alexander E. Lipka (2016)


**K_Chr_Generalized.R** contains a "K_Chr function" created by Alex Lipka and Angela Chen. It creates "K_Chr kinship matrix" and "K_Chr principal components of population structure" and then runs GWAS. Details of this method can be found in our manuscript. This function can be adapted to analyze any datasets.  

It reads in genotypic data in hapmap format, with a separate file for each chromosome, and phenotypic data. Please see the GAPIT user manual for a description of how to format the input file. Manual here: http://www.zzlab.net/GAPIT/gapit_help_document.pdf

It generates GWAS results in a table, containing important information such as p values, FDR adjusted p values, etc. It also produces Manhattan plots that show the corresponding genomic position of a significant SNP identified.

**Run_K_Chr_Generalized.R** is an execution file that runs the K_Chr function. Please note that you will need to read in phenotypic data here. 

The folder **GAPIT_Demo_Data** contains genotypic data and phenotypic data as an example for running **Run_K_Chr_Generalized.R**.
