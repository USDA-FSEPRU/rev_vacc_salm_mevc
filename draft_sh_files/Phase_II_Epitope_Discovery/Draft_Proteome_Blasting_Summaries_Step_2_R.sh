##Preping R environment

#Install packages not in default R in SCINet
#install.packages(c("tidyverse")

#Load libraries
library(tidyverse)

#Load bash variables into R environment
args <- commandArgs()

#Load the bash variable "epitope_unique_id" to create an R variable
Unique_Id <- args[7]

#Create a variable for the path from the bash variable "summary_file" 
path <- paste0("Infantis_splitting/", args[6])

#Upload results from summarizing blasting proteome results for one epitope
proteomes_epitope_homology <- read.delim(path, header=FALSE, sep=" ")

#Determine percentage of proteomes with 100% homology
percent_positive_epitope_homology <- sum(proteomes_epitope_homology$V3=="TRUE")/14052*100   #Calculate the percentage of TRUEs across all proteomes

#Create a series of TRUE/FALSE statements regarding how high this percentage is
if (percent_positive_epitope_homology == 100){                                                                  #If the percentage above is = 100...
    Homology = "TRUE"                                                                                             #Say Homology is TRUE
  } else {                                                                                                        #Otherwise...
    Homology = "FALSE"                                                                                            #Say homology is FALSE
  }
  
  if (percent_positive_epitope_homology > 99){                                                                  #If the percentage above is > 99...
    Partial_99_Homology = "TRUE"                                                                                             #Say Homology is TRUE
  } else {                                                                                                        #Otherwise...
    Partial_99_Homology = "FALSE"                                                                                            #Say homology is FALSE
  }  
  
  if (percent_positive_epitope_homology > 98){                                                                  #If the percentage above is > 98...
    Partial_98_Homology = "TRUE"                                                                                             #Say Homology is TRUE
  } else {                                                                                                        #Otherwise...
    Partial_98_Homology = "FALSE"                                                                                            #Say homology is FALSE
  }  
  
  if (percent_positive_epitope_homology > 97){                                                                  #If the percentage above is > 97...
    Partial_97_Homology = "TRUE"                                                                                             #Say Homology is TRUE
  } else {                                                                                                        #Otherwise...
    Partial_97_Homology = "FALSE"                                                                                            #Say homology is FALSE
  }   
  
  if (percent_positive_epitope_homology > 96){                                                                  #If the percentage above is > 96...
    Partial_96_Homology = "TRUE"                                                                                             #Say Homology is TRUE
  } else {                                                                                                        #Otherwise...
    Partial_96_Homology = "FALSE"                                                                                            #Say homology is FALSE
  }  
  
  if (percent_positive_epitope_homology > 95){                                                                  #If the percentage above is > 95...
    Partial_95_Homology = "TRUE"                                                                                             #Say Homology is TRUE
  } else {                                                                                                        #Otherwise...
    Partial_95_Homology = "FALSE"                                                                                            #Say homology is FALSE
  }  
  
  if (percent_positive_epitope_homology > 90){                                                                  #If the percentage above > = 90...
    Partial_90_Homology = "TRUE"                                                                                             #Say Homology is TRUE
  } else {                                                                                                        #Otherwise...
    Partial_90_Homology = "FALSE"                                                                                            #Say homology is FALSE
  }  


#Put all the information into a dataframe
proteomes_epitope_homology_summary <- data.frame(Unique_Id, percent_positive_epitope_homology, Homology, Partial_99_Homology, Partial_98_Homology, Partial_97_Homology, Partial_96_Homology, Partial_95_Homology, Partial_90_Homology)  #Save all elements above in the summary dataframe


#Add column names
colnames(proteomes_epitope_homology_summary) <- c("Unique_Id", "Epitope_vs_Infantis_proteomes_Homology_Percentage", "Epitope_vs_Infantis_proteomes_100_percent_Homology", "Epitope_vs_Infantis_proteomes_99_percent_Homology", "Epitope_vs_Infantis_proteomes_98_percent_Homology", "Epitope_vs_Infantis_proteomes_97_percent_Homology", "Epitope_vs_Infantis_proteomes_96_percent_Homology", "Epitope_vs_Infantis_proteomes_95_percent_Homology", "Epitope_vs_Infantis_proteomes_90_percent_Homology")

#Split Unique_Id column into Identity and Extra bits
proteomes_epitope_homology_summary <- separate(data = proteomes_epitope_homology_summary, col = Unique_Id, into = c("Identity", "Extra"), sep = "_", extra = "merge", remove = FALSE)

#Write out the summary for this epitope
write_delim(proteomes_epitope_homology_summary, paste0("Infantis_faas/epitope_blasting_summaries_Step_2/", Unique_Id, "_proteomes_epitope_homology_summary.txt"))
