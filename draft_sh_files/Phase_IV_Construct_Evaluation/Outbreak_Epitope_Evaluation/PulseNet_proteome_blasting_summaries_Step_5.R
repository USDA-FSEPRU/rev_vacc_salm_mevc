##Preping R environment

#Install packages not in default R in SCINet
#install.packages(c("tidyverse")

#Load libraries
library(tidyverse)
library(foreach)

#Load bash variables into R environment
args <- commandArgs()

#Load the bash variable "epitope_unique_id" to create an R variable
Unique_Id <- args[7]

#Create a variable for the path from the bash variable "summary_file" 
path <- paste0("PulseNet_splitting/", args[6])

#Upload results from summarizing blasting proteome results for one epitope
proteomes_epitope_homology <- read.delim(path, header=FALSE, sep=" ")

#Upload the PulseNet SRR Accession Number to Isolation Source dataframe
PulseNet_All_summ_CheckM_Pass_Iso_Src <- read.delim("PulseNet_All_summ_CheckM_Pass_Iso_Src.tsv")

#Combine dataframes
proteomes_epitope_homology <-
  merge(proteomes_epitope_homology, 
        PulseNet_All_summ_CheckM_Pass_Iso_Src,
        by.x="V2", 
        by.y="SRR_id")

#Check for number of complete homology hits for each Isolation Source
proteomes_epitope_homology_isosrc <- foreach(IsoSrc=unique(proteomes_epitope_homology$Summarized.isolation.source), .combine = rbind) %do% {
  subset_isosrc = filter(proteomes_epitope_homology, Summarized.isolation.source == IsoSrc)
  full_0_epitope_percentage <- sum(subset_isosrc$V6=="TRUE")/nrow(subset_isosrc)*100
  data.frame(IsoSrc, full_0_epitope_percentage) #Bind all elements above to data frame
}

#Change the column names
colnames(proteomes_epitope_homology_isosrc) <- c("IsoSrc", Unique_Id)

#Make IsoSrc the rownames
proteomes_epitope_homology_isosrc <- column_to_rownames(proteomes_epitope_homology_isosrc,"IsoSrc")

#Transpose the dataframe
proteomes_epitope_homology_isosrc <- data.frame(t(proteomes_epitope_homology_isosrc))

#Put IsoSrc columns in alphabetical order
proteomes_epitope_homology_isosrc <- proteomes_epitope_homology_isosrc %>% select(order(colnames(proteomes_epitope_homology_isosrc)))

#Make rowname into Unique_Id column
proteomes_epitope_homology_isosrc <- rownames_to_column(proteomes_epitope_homology_isosrc,"Unique_Id")

#Write out the summary for this epitope
write_tsv(proteomes_epitope_homology_isosrc, paste0("epitope_blasting_results_summaries_Step_5/", Unique_Id, "_proteomes_epitope_homology_isosrc_summary.tsv"))

