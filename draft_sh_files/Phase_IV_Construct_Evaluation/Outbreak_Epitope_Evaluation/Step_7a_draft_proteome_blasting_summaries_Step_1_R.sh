library(dplyr)
library(foreach)
library(readr)

args <- commandArgs()

#Load the bash variable "proteome_name" to create an R variable
Proteome_Name <- args[7]

#Create a variable for the path from the bash variable "blast_file" 
path <- paste0("epitope_blasting_results/", args[6])

#Load up the blast results from one proteome
blast_results <- read.delim(path, header=FALSE)

#Make sure pident is a numerical class
blast_results$V3 <- as.numeric(blast_results$V3)

#Make sure qcovhsp is a numerical class
blast_results$V15 <- as.numeric(blast_results$V15)

#Make sure length is a numerical class
blast_results$V4 <- as.numeric(blast_results$V4)

#Make a for loop to determine if a epitope had 100% identity to the proteome or not
blast_results_epitope_homology <- foreach(Unique_Id=unique(blast_results$V1), .combine = rbind) %do% {
  subset_epitope = filter(blast_results, V1 == Unique_Id)
  pident_homology = ifelse(max(subset_epitope$V3) == 100, TRUE, FALSE)
  cov_homology = ifelse(max(subset_epitope$V15) == 100, TRUE, FALSE)
  max_length = max(subset_epitope$V4)
  subset_epitope$full_homology = ifelse(subset_epitope$V3 == 100 & subset_epitope$V15 == 100, "Yes", "No")
  full_0_homology = ifelse("Yes" %in% subset_epitope$full_homology, TRUE, FALSE)
  subset_epitope$partial_homology = ifelse(subset_epitope$V3 >= 88 & subset_epitope$V15 >= 88, "Yes", "No")
  partial_1_homology = ifelse("Yes" %in% subset_epitope$partial_homology, TRUE, FALSE)
  max_pident = max(subset_epitope$V3)
  max_qcovhsp = max(subset_epitope$V15)
  data.frame(Unique_Id, Proteome_Name, pident_homology, max_length, cov_homology, full_0_homology, partial_1_homology, max_pident, max_qcovhsp) #Bind all elements above to data frame
}

#Write out the resulting dataframe
write_delim(blast_results_epitope_homology, paste0("epitope_blasting_results_summaries_Step_1/", Proteome_Name, "_proteome_epitope_homology.txt"))

