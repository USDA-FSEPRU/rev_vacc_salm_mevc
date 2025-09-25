##Preping R environment

#Install packages not in default R in SCINet
#install.packages(c("tidyverse", "doParallel"))

#Load Libraries
library(tidyverse)
library(doParallel)
library(parallel)
library(foreach)

##Hadar

#Upload results from blasting the proteomes
Hadar_proteomes <- read.delim("Hadar_faas/Hadar_proteomes_epitope_blasting_results_070723_R.txt", header = FALSE)

#Create new columns separating out the Proteome_Name and Unique_Id
Hadar_proteomes <- separate(data = Hadar_proteomes, col = V1, into = c("Proteome_Name", "Unique_Id"), sep = ";", remove = FALSE)

#Change column names to something readable
colnames(Hadar_proteomes) <- c("Proteome_Name_Unique_Id", "Proteome_Name", "Unique_Id", "pident")

#Make sure pident is a numerical class
Hadar_proteomes$pident <- as.numeric(Hadar_proteomes$pident)

#Prep for running in parallel

#Detect number of available cores and create cluster
cl <- parallel::makeCluster(24)

# Activate cluster for foreach library
doParallel::registerDoParallel(cl)

#Make a for loop to determine if a epitope had 100% identity to a strain or not
Hadar_proteomes_epitope_homology <- foreach(proteome_epitope=unique(Hadar_proteomes$Proteome_Name_Unique_Id), .combine = rbind, .packages=c('tidyverse')) %dopar% {
  subset_epitope = filter(Hadar_proteomes, Proteome_Name_Unique_Id == proteome_epitope)
  Proteome_Name <- subset_epitope[1, "Proteome_Name"]                                                                    #Save the Proteome_Name
  Unique_Id <- subset_epitope[1, "Unique_Id"]
  if (max(subset_epitope$pident) == 100){                                                                            #If the max percent idenity is = 100...
    homology = "TRUE"                                                                                                #Say Homology is TRUE...
  } else {                                                                                                           #Otherwise...
    homology = "FALSE"                                                                                               #Say Homology is FALSE
  }
  data.frame(Unique_Id, Proteome_Name, homology) #Bind three elements above to data frame
  
}

#Change column names to something readable
colnames(Hadar_proteomes_epitope_homology) <- c("Unique_Id", "Proteome_Name", "Homology")

#Create a blank dataframe to store a per epitope summary of percent identity to all proteomes
Hadar_proteomes_epitope_homology_summary <- data.frame()

#Create a loop to determine what percentage of proteomes the epitope had 100% identity with
Hadar_proteomes_epitope_homology_summary <- foreach(uniqueid = unique(Hadar_proteomes_epitope_homology$Unique_Id), .combine = rbind, .packages=c('tidyverse')) %dopar%{                                           #For each Unique_Id
  subset_epitope = filter(Hadar_proteomes_epitope_homology, Unique_Id == uniqueid)                               #Subset the dataframe to just that Unique_Id
  percent_positive_epitope_homology <- sum(subset_epitope$Homology=="TRUE")/length(subset_epitope$Homology)*100   #Calculate the percentage of TRUEs across all proteomes
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
   data.frame(uniqueid, percent_positive_epitope_homology, Homology, Partial_99_Homology, Partial_98_Homology, Partial_97_Homology, Partial_96_Homology, Partial_95_Homology, Partial_90_Homology)  #Save all elements above in the summary dataframe
}

#Change column names to something readable
colnames(Hadar_proteomes_epitope_homology_summary) <- c("Unique_Id", "Epitope_vs_Hadar_proteomes_Homology_Percentage", "Epitope_vs_Hadar_proteomes_100_percent_Homology", "Epitope_vs_Hadar_proteomes_99_percent_Homology", "Epitope_vs_Hadar_proteomes_98_percent_Homology", "Epitope_vs_Hadar_proteomes_97_percent_Homology", "Epitope_vs_Hadar_proteomes_96_percent_Homology", "Epitope_vs_Hadar_proteomes_95_percent_Homology", "Epitope_vs_Hadar_proteomes_90_percent_Homology")

#Split Unique_Id column into Identity and Extra bits
Hadar_proteomes_epitope_homology_summary <- separate(data = Hadar_proteomes_epitope_homology_summary, col = Unique_Id, into = c("Identity", "Extra"), sep = "_", extra = "merge", remove = FALSE)

#Check distribution of TRUEs and FALSEs
Hadar_proteomes_epitope_homology_counts_summary <- as.data.frame(cbind(t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_100_percent_Homology)), t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_99_percent_Homology)), t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_95_percent_Homology)), t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_90_percent_Homology))))

write_delim(Hadar_proteomes_epitope_homology, "Hadar_proteomes_epitope_homology.txt")
write_delim(Hadar_proteomes_epitope_homology_summary, "Hadar_proteomes_epitope_homology_summary_v2.txt")
write_delim(Hadar_proteomes_epitope_homology_counts_summary, "Hadar_proteomes_epitope_homology_counts_summary.txt")
save.image("Hadar_V2.RData")
