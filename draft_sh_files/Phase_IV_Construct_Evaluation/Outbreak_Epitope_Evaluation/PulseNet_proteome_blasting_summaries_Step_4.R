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

#Upload the PulseNet SRR Accession Number to Serovar dataframe
PulseNet_All_summ_CheckM_Pass_Serovar <- read.delim("PulseNet_All_summ_CheckM_Pass_Serovar.tsv")

#Combine dataframes
proteomes_epitope_homology <-
  merge(proteomes_epitope_homology, 
        PulseNet_All_summ_CheckM_Pass_Serovar,
        by.x="V2", 
        by.y="SRR_id")

#Change NAs to None in Serotype_wgs column
proteomes_epitope_homology$Serotype_wgs[is.na(proteomes_epitope_homology$Serotype_wgs)] = "None"

#Check for number of complete homology hits for each serovar
proteomes_epitope_homology_serovar <- foreach(Serovar=unique(proteomes_epitope_homology$Serotype_wgs), .combine = rbind) %do% {
  subset_serovar = filter(proteomes_epitope_homology, Serotype_wgs == Serovar)
  full_0_epitope_count<- sum(subset_serovar$V6=="TRUE")
  Unique_Id = subset_serovar$V1
  data.frame(Serovar, full_0_epitope_count) #Bind all elements above to data frame
}

#Change the column names
colnames(proteomes_epitope_homology_serovar) <- c("Serovar", Unique_Id)

#Make serovar the rownames
proteomes_epitope_homology_serovar <- column_to_rownames(proteomes_epitope_homology_serovar,"Serovar")

#Transpose the dataframe
proteomes_epitope_homology_serovar <- data.frame(t(proteomes_epitope_homology_serovar))

#Put Serovar columns in alphabetical order
proteomes_epitope_homology_serovar <- proteomes_epitope_homology_serovar %>% select(order(colnames(proteomes_epitope_homology_serovar)))

#Make rowname into Unique_Id column
proteomes_epitope_homology_serovar <- rownames_to_column(proteomes_epitope_homology_serovar,"Unique_Id")

#Write out the summary for this epitope
write_tsv(proteomes_epitope_homology_serovar, paste0("epitope_blasting_results_summaries_Step_4/", Unique_Id, "_proteomes_epitope_homology_serovar_counts.tsv"))

