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
PulseNet_All_summ_CheckM_Pass_Otbk_Src <- read.delim("PulseNet_All_summ_CheckM_Pass_Otbk_Src.tsv")

#Combine dataframes
proteomes_epitope_homology <-
  merge(proteomes_epitope_homology, 
        PulseNet_All_summ_CheckM_Pass_Otbk_Src,
        by.x="V2", 
        by.y="SRR_id")

#Check for number of complete homology hits for each otbksrc
proteomes_epitope_homology_otbksrc <- foreach(OtbkSrc=unique(proteomes_epitope_homology$Summarized.outbreak.source), .combine = rbind) %do% {
  subset_otbksrc = filter(proteomes_epitope_homology, Summarized.outbreak.source == OtbkSrc)
  full_0_epitope_count<- sum(subset_otbksrc$V6=="TRUE")
  data.frame(OtbkSrc, full_0_epitope_count) #Bind all elements above to data frame
}

#Change the column names
colnames(proteomes_epitope_homology_otbksrc) <- c("OtbkSrc", Unique_Id)

#Make OtbkSrc the rownames
proteomes_epitope_homology_otbksrc <- column_to_rownames(proteomes_epitope_homology_otbksrc,"OtbkSrc")

#Transpose the dataframe
proteomes_epitope_homology_otbksrc <- data.frame(t(proteomes_epitope_homology_otbksrc))

#Put OtbkSrc columns in alphabetical order
proteomes_epitope_homology_otbksrc <- proteomes_epitope_homology_otbksrc %>% select(order(colnames(proteomes_epitope_homology_otbksrc)))

#Make rowname into Unique_Id column
proteomes_epitope_homology_otbksrc <- rownames_to_column(proteomes_epitope_homology_otbksrc,"Unique_Id")

#Write out the summary for this epitope
write_tsv(proteomes_epitope_homology_otbksrc, paste0("epitope_blasting_results_summaries_Step_8/", Unique_Id, "_proteomes_epitope_homology_otbksrc_counts.tsv"))

