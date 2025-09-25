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
path <- paste0("PulseNet_splitting/", args[6])

#Upload results from summarizing blasting proteome results for one epitope
proteomes_epitope_homology <- read.delim(path, header=FALSE, sep=" ")

##100% Percent Identity with Varying Coverage

#Determine percentage of proteomes with 100% pident
pident_epitope_percentage <- sum(proteomes_epitope_homology$V3=="TRUE")/4*100   #Calculate the percentage of TRUEs across all proteomes

#Create a series of TRUE/FALSE statements regarding how high this percentage is
pident_ET_100 <- ifelse(pident_epitope_percentage == 100, TRUE, FALSE) #If the percentage is equal to 100, say homology is TRUE, otherwise say Homology is FALSE
pident_GT_99 <- ifelse(pident_epitope_percentage > 99, TRUE, FALSE) #If the percentage is greater than to 99, say homology is TRUE, otherwise say Homology is FALSE
pident_GT_98 <- ifelse(pident_epitope_percentage > 98, TRUE, FALSE) #If the percentage is greater than to 98, say homology is TRUE, otherwise say Homology is FALSE
pident_GT_97 <- ifelse(pident_epitope_percentage > 97, TRUE, FALSE) #If the percentage is greater than to 97, say homology is TRUE, otherwise say Homology is FALSE
pident_GT_96 <- ifelse(pident_epitope_percentage > 96, TRUE, FALSE) #If the percentage is greater than to 96, say homology is TRUE, otherwise say Homology is FALSE
pident_GT_95 <- ifelse(pident_epitope_percentage > 95, TRUE, FALSE) #If the percentage is greater than to 95, say homology is TRUE, otherwise say Homology is FALSE
pident_GT_90 <- ifelse(pident_epitope_percentage > 90, TRUE, FALSE) #If the percentage is greater than to 90, say homology is TRUE, otherwise say Homology is FALSE

##100% Coverage with Varying Percent Identity

#Determine percentage of proteomes with 100% coverage
qcovhsp_epitope_percentage <- sum(proteomes_epitope_homology$V5=="TRUE")/4*100   #Calculate the percentage of TRUEs across all proteomes

#Create a series of TRUE/FALSE statements regarding how high this percentage is
qcovhsp_ET_100 <- ifelse(qcovhsp_epitope_percentage == 100, TRUE, FALSE) #If the percentage is equal to 100, say homology is TRUE, otherwise say Homology is FALSE
qcovhsp_GT_99 <- ifelse(qcovhsp_epitope_percentage > 99, TRUE, FALSE) #If the percentage is greater than to 99, say homology is TRUE, otherwise say Homology is FALSE
qcovhsp_GT_98 <- ifelse(qcovhsp_epitope_percentage > 98, TRUE, FALSE) #If the percentage is greater than to 98, say homology is TRUE, otherwise say Homology is FALSE
qcovhsp_GT_97 <- ifelse(qcovhsp_epitope_percentage > 97, TRUE, FALSE) #If the percentage is greater than to 97, say homology is TRUE, otherwise say Homology is FALSE
qcovhsp_GT_96 <- ifelse(qcovhsp_epitope_percentage > 96, TRUE, FALSE) #If the percentage is greater than to 96, say homology is TRUE, otherwise say Homology is FALSE
qcovhsp_GT_95 <- ifelse(qcovhsp_epitope_percentage > 95, TRUE, FALSE) #If the percentage is greater than to 95, say homology is TRUE, otherwise say Homology is FALSE
qcovhsp_GT_90 <- ifelse(qcovhsp_epitope_percentage > 90, TRUE, FALSE) #If the percentage is greater than to 90, say homology is TRUE, otherwise say Homology is FALSE

##100% Coverage and Percent Identity

#Determine percentage of proteomes with 100% coverage and percent identity
full_0_epitope_percentage <- sum(proteomes_epitope_homology$V6=="TRUE")/4*100   #Calculate the percentage of TRUEs across all proteomes

#Create a series of TRUE/FALSE statements regarding how high this percentage is
full_0_ET_100 <- ifelse(full_0_epitope_percentage == 100, TRUE, FALSE) #If the percentage is equal to 100, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_99 <- ifelse(full_0_epitope_percentage > 99, TRUE, FALSE) #If the percentage is greater than to 99, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_98 <- ifelse(full_0_epitope_percentage > 98, TRUE, FALSE) #If the percentage is greater than to 98, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_97 <- ifelse(full_0_epitope_percentage > 97, TRUE, FALSE) #If the percentage is greater than to 97, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_96 <- ifelse(full_0_epitope_percentage > 96, TRUE, FALSE) #If the percentage is greater than to 96, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_95 <- ifelse(full_0_epitope_percentage > 95, TRUE, FALSE) #If the percentage is greater than to 95, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_90 <- ifelse(full_0_epitope_percentage > 90, TRUE, FALSE) #If the percentage is greater than to 90, say homology is TRUE, otherwise say Homology is FALSE

#Determine percentage of proteomes with >88% coverage and percent identity
partial_1_epitope_percentage <- sum(proteomes_epitope_homology$V7=="TRUE")/4*100   #Calculate the percentage of TRUEs across all proteomes

#Create a series of TRUE/FALSE statements regarding how high this percentage is
partial_1_ET_100 <- ifelse(partial_1_epitope_percentage == 100, TRUE, FALSE) #If the percentage is equal to 100, say homology is TRUE, otherwise say Homology is FALSE
partial_1_GT_99 <- ifelse(partial_1_epitope_percentage > 99, TRUE, FALSE) #If the percentage is greater than to 99, say homology is TRUE, otherwise say Homology is FALSE
partial_1_GT_98 <- ifelse(partial_1_epitope_percentage > 98, TRUE, FALSE) #If the percentage is greater than to 98, say homology is TRUE, otherwise say Homology is FALSE
partial_1_GT_97 <- ifelse(partial_1_epitope_percentage > 97, TRUE, FALSE) #If the percentage is greater than to 97, say homology is TRUE, otherwise say Homology is FALSE
partial_1_GT_96 <- ifelse(partial_1_epitope_percentage > 96, TRUE, FALSE) #If the percentage is greater than to 96, say homology is TRUE, otherwise say Homology is FALSE
partial_1_GT_95 <- ifelse(partial_1_epitope_percentage > 95, TRUE, FALSE) #If the percentage is greater than to 95, say homology is TRUE, otherwise say Homology is FALSE
partial_1_GT_90 <- ifelse(partial_1_epitope_percentage > 90, TRUE, FALSE) #If the percentage is greater than to 90, say homology is TRUE, otherwise say Homology is FALSE


#Put all the information into a dataframe
proteomes_epitope_homology_summary <- data.frame(Unique_Id, pident_epitope_percentage, pident_ET_100, pident_GT_99, pident_GT_98, pident_GT_97, pident_GT_96, pident_GT_95, pident_GT_90,
                                                 qcovhsp_epitope_percentage, qcovhsp_ET_100, qcovhsp_GT_99, qcovhsp_GT_98, qcovhsp_GT_97, qcovhsp_GT_96, qcovhsp_GT_95, qcovhsp_GT_90,
                                                 full_0_epitope_percentage, full_0_ET_100, full_0_GT_99, full_0_GT_98, full_0_GT_97, full_0_GT_96, full_0_GT_95, full_0_GT_90,
                                                 partial_1_epitope_percentage, partial_1_ET_100, partial_1_GT_99, partial_1_GT_98, partial_1_GT_97, partial_1_GT_96, partial_1_GT_95, partial_1_GT_90)  #Save all elements above in the summary dataframe

#Add column names
colnames(proteomes_epitope_homology_summary) <- c("Unique_Id", "Epitope_vs_PN_Pident_Percentage", "Pident_ET_100", "Pident_GT_99", "Pident_GT_98", "Pident_GT_97", "Pident_GT_96", "Pident_GT_95", "Pident_GT_90",
                                                  "Epitope_vs_PN_Qcovhsp_Percentage", "Qcovhsp_ET_100", "Qcovhsp_GT_99", "Qcovhsp_GT_98", "Qcovhsp_GT_97", "Qcovhsp_GT_96", "Qcovhsp_GT_95", "Qcovhsp_GT_90",
                                                  "Epitope_vs_PN_Full_0_Percentage", "Full_0_ET_100", "Full_0_GT_99", "Full_0_GT_98", "Full_0_GT_97", "Full_0_GT_96", "Full_0_GT_95", "Full_0_GT_90",
                                                  "Epitope_vs_PN_Partial_1_Percentage", "Partial_1_ET_100", "Partial_1_GT_99", "Partial_1_GT_98", "Partial_1_GT_97", "Partial_1_GT_96", "Partial_1_GT_95", "Partial_1_GT_90")

#Split Unique_Id column into Identity and Extra bits
proteomes_epitope_homology_summary <- separate(data = proteomes_epitope_homology_summary, col = Unique_Id, into = c("Identity", "Extra"), sep = "_", extra = "merge", remove = FALSE)

#Write out the summary for this epitope
write_delim(proteomes_epitope_homology_summary, paste0("epitope_blasting_results_summaries_Step_2/", Unique_Id, "_proteomes_epitope_homology_summary.txt"))
