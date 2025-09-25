library(dplyr)
library(foreach)
library(readr)

args <- commandArgs()

#Load the bash variable "proteome_name" to create an R variable
Proteome_Name <- args[7]

#Create a variable for the path from the bash variable "blast_file" 
path <- paste0("Enteritidis_faas/epitope_blasting_results/", args[6])

#Load up the blast results from one proteome
blast_results <- read.delim(path, header=FALSE)

#Make sure pident is a numerical class
blast_results$V3 <- as.numeric(blast_results$V3)

#Make a for loop to determine if a epitope had 100% identity to the proteome or not
blast_results_epitope_homology <- foreach(Unique_Id=unique(blast_results$V1), .combine = rbind) %do% {
  subset_epitope = filter(blast_results, V1 == Unique_Id)
  if (max(subset_epitope$V3) == 100){                                                                            #If the max percent idenity is = 100...
    homology = "TRUE"                                                                                                #Say Homology is TRUE...
  } else {                                                                                                           #Otherwise...
    homology = "FALSE"                                                                                               #Say Homology is FALSE
  }
  data.frame(Unique_Id, Proteome_Name, homology) #Bind three elements above to data frame
  
}

#Write out the resulting dataframe
write_delim(blast_results_epitope_homology, paste0("Enteritidis_faas/epitope_blasting_summaries_Step_1/", Proteome_Name, "_proteome_epitope_homology.txt"))
