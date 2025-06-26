#Set overall working directory
setwd("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Projects_Reverse_Vaccinology/UK1/cleaned_up")
#Load libraries
library(tidyverse)
library(Biostrings)
library(ape)
library(seqinr)


#See Reverse_Vaccinology_UK1_Scripts.txt for initial information regarding getting proteome

#See you back here once you are finished with Vaxign2!

####PHASE I SUBTRACTIVE PROTEOMICS####

#...
#Importing information from Excel/Notepad++...
##Number of Proteins in Original UK1 GenBank Proteome = 4555###
#...

#Set working directory to Phase specific folder
setwd("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Projects_Reverse_Vaccinology/UK1/cleaned_up/Phase_I_Subtractive_Proteomics")
#Upload the results as a tab deliminated file
Full_UK1_Vaxign2 <- read.delim("Full_UK1_Vaxign2_Results.txt") #4555 x 7 dataframe

#Make a vector of the localizations you want to keep
#These localizations were used in other papers
#Although most papers just focus on Outer Membrane Proteins, we wanted to have a MESV with multiple localizations
Localizations_to_Keep <- c("Extracellular", "Outer Membrane", "Periplasmic")

#The rest of the filters (transmembrance helices <=1, adhesion probability >=0.51, and Vaxign ML score >=90 are defaults from Vaxign2's Precomputed Query server and documentation)
Filtered_UK1_Vaxign2 <- filter(Full_UK1_Vaxign2, Trans.membrane.Helices <= 1 & Adhesin.Probability >= 0.51 & Vaxign.ML.Score >= 90 & Localization %in% Localizations_to_Keep)

#Copy table to Excel to track changes
clipr::write_clip(Filtered_UK1_Vaxign2)

##Number of Proteins post-Vaxign2 filtering = 190###

####Get Additional Annotations from Similar Salmonella Proteins####

#See the Reverse_Vaccinology_UK1_Scripts.txt for further instructions, search for header above "###Get Additional Annotations from Similar Salmonella Proteins"
#Good luck with blasting!
#...
#Welcome back!

UK1_ReAnno <- read.delim("UK1_ReAnno_tabular_060923_edited.txt", header = TRUE)

#Get an idea of what the minimum percent identities were for the proteins
#Good way to sanity check that you have right number of proteins too 
#190 for this study
UK1_ReAnno_pident_mins = data.frame()
for (protein in unique(UK1_ReAnno$qaccver)) {
  subset_protein = filter(UK1_ReAnno, qaccver == protein)
  pindent_min <- min(subset_protein$pident)
  UK1_ReAnno_pident_mins <- rbind(UK1_ReAnno_pident_mins, c(protein, pindent_min))
}

#Remove intermediates
remove(subset_protein)
remove(protein)
remove(pindent_min)

#Filter any percent identities lower than 98%, want a high level of homology
#Per https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10143441/
#Means we are keeping either the first 500 hits (BLAST default) or those with a percent identity >98%
UK1_ReAnno <- filter(UK1_ReAnno, pident >= 98)

#Gather information from each row and summarize it into a useful table
#This'll take some time
UK1_ReAnno_summary = data.frame()
for(row in 1:nrow(UK1_ReAnno)){
  protein_query <- UK1_ReAnno[row, "qaccver"]         # Save the protein query name
  subject_description <- UK1_ReAnno[row, "stitle"]    #save the subject description
  gene_name <- word(subject_description, -2)          #Extract out the last word at the end of the description
  if (str_detect(gene_name,"[[:upper:]]") == TRUE){   #If it has an upper case letter then go into loop
    if (str_detect(gene_name, "[0-9]") == FALSE) {    #If it does not have a number...
      gene_name = gene_name                           #Save the word as the assumed gene name
    }
    else {                                            #If it does have a number...
      gene_name = NA                                  #Put NA
    }
  } else {                                            #If it does not have an upper case letter...
    gene_name = NA                                    #Put NA
  }
  UK1_ReAnno_summary <- rbind(UK1_ReAnno_summary, c(protein_query, subject_description, gene_name)) #Save the three elements per row as a dataframe
}

#Remove intermediates
remove(row)
remove(protein_query)
remove(subject_description)
remove(gene_name)

#Change the column names to something useful
colnames(UK1_ReAnno_summary) <- c("GenBank.Accession", "Subject.Description", "Subject.Gene.Name")

#Summarize the above information on a per protein basis
UK1_ReAnno_gene_name_summary = data.frame()
for (protein in unique(UK1_ReAnno_summary$GenBank.Accession)) {                   #For each unique protein...
  subset_protein = filter(UK1_ReAnno_summary, GenBank.Accession == protein)       #Subset the main dataframe to just that protein
  gene_names <- unique(subset_protein$Subject.Gene.Name)                      #Keep a list of the unique subject gene names
  gene_names <- gene_names[!is.na(gene_names)]                                #Get rid of any NAs in the list
  gene_descriptions <- unique(subset_protein$Subject.Description)             #Keep a list of the unique subject descriptions
  UK1_ReAnno_gene_name_summary <- rbind(UK1_ReAnno_gene_name_summary, c(protein, toString(gene_descriptions), toString(gene_names))) #Save the three elements as a dataframe, toString makes the list into a string to put into a column
}

#Remove intermediates
remove(subset_protein)
remove(gene_names)
remove(gene_descriptions)

#Change colnames to something useful
colnames(UK1_ReAnno_gene_name_summary) <- c("GenBank.Accession", "Similar.Salmonella.Proteins.Descriptions", "Similar.Salmonella.Proteins.Protein.Names")

#Change any NAs to " -"
UK1_ReAnno_gene_name_summary[is.na(UK1_ReAnno_gene_name_summary)] = " -"

#Merge with Vaxign2 results by GenBank.Accession
Filtered_UK1_Vaxign2_ReAnno <- merge(Filtered_UK1_Vaxign2, UK1_ReAnno_gene_name_summary)

#Copy down to Excel for manipulation
clipr::write_clip(Filtered_UK1_Vaxign2_ReAnno)

#Visually check and edit the resulting table
#Keep in mind that R code assumes no there are no digits [0-9] in the protein name as well as at least one capital letter

####Removal of Plasmid, Flagellar, and LPS Proteins####

#See Reverse_Vaccinology_UK1_Scripts.txt for initial information regarding next filtering steps

#See ya after you are done with Protein Blasting!
#...

#Importing information from Excel/Notepad++...
###Non-flagellar, -LPS, and -plasmid UK1 Proteins = 167###
###Antigenic Proteins = 127###

#...
#Welcome back!

####Negative homology####

#Set working directory to location of all the protein blasting data
setwd("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Projects_Reverse_Vaccinology/UK1/cleaned_up/Phase_I_Subtractive_Proteomics/Protein_Blasting")
##Gallus gallus
Gallus_gallus_protein_blasting_table<- read.table("UK1_Gallus_gallus_tabular_070323.txt")        #Upload your tabular results
Gallus_gallus_protein_homology = data.frame()                                                #Create a blank dataframe to store summarized results
for (protein in unique(Gallus_gallus_protein_blasting_table$V1)) {                                    #For every unique protein query (UK1 Accession)...
  subset_protein = filter(Gallus_gallus_protein_blasting_table, V1 == protein)                        #Create a subset of the original dataframe with just that protein query
  if (min(subset_protein$V11) < 1e-4){                                                   #If the minimum evalue is < 1e-4 then...
    homology = "True"                                                                    #State "True" aka protein is homologous to host proteins
  } else {
    homology = "False"                                                                   #Otherwise state "False" aka protein is nonhomologous to host proteins
  }
  Gallus_gallus_protein_homology <- rbind(Gallus_gallus_protein_homology, c(protein, homology))  #Save protein query and summary as a line in summarized dataframe
}
colnames(Gallus_gallus_protein_homology) <- c("GenBank.Accession", "Gallus.gallus.Homology") #Change column names to something readable

##Meleagris gallopavo
Meleagris_gallopavo_protein_blasting_table<- read.table("UK1_Meleagris_gallopavo_tabular_070323.txt")
Meleagris_gallopavo_protein_homology = data.frame()
for (protein in unique(Meleagris_gallopavo_protein_blasting_table$V1)) {
  subset_protein = filter(Meleagris_gallopavo_protein_blasting_table, V1 == protein)
  if (min(subset_protein$V11) < 1e-4){
    homology = "True"
  } else {
    homology = "False"
  }
  Meleagris_gallopavo_protein_homology <- rbind(Meleagris_gallopavo_protein_homology, c(protein, homology))
}
colnames(Meleagris_gallopavo_protein_homology) <- c("GenBank.Accession", "Meleagris.gallopavo.Homology")

##Sus scrofa
Sus_scrofa_protein_blasting_table <- read.table("UK1_Sus_scrofa_tabular_070323.txt")
Sus_scrofa_protein_homology = data.frame()
for (protein in unique(Sus_scrofa_protein_blasting_table$V1)) {
  subset_protein = filter(Sus_scrofa_protein_blasting_table, V1 == protein)
  if (min(subset_protein$V11) < 1e-4){
    homology = "True"
  } else {
    homology = "False"
  }
  Sus_scrofa_protein_homology <- rbind(Sus_scrofa_protein_homology, c(protein, homology))
}
colnames(Sus_scrofa_protein_homology) <- c("GenBank.Accession", "Sus.scrofa.Homology")

##Mus musculus
Mus_musculus_protein_blasting_table <- read.table("UK1_Mus_musculus_tabular_070323.txt")
Mus_musculus_protein_homology = data.frame()
for (protein in unique(Mus_musculus_protein_blasting_table$V1)) {
  subset_protein = filter(Mus_musculus_protein_blasting_table, V1 == protein)
  if (min(subset_protein$V11) < 1e-4){
    homology = "True"
  } else {
    homology = "False"
  }
  Mus_musculus_protein_homology <- rbind(Mus_musculus_protein_homology, c(protein, homology))
}
colnames(Mus_musculus_protein_homology) <- c("GenBank.Accession", "Mus.musculus.Homology")

##Bos taurus
Bos_taurus_protein_blasting_table <- read.table("UK1_Bos_taurus_tabular_070323.txt")
Bos_taurus_protein_homology = data.frame()
for (protein in unique(Bos_taurus_protein_blasting_table$V1)) {
  subset_protein = filter(Bos_taurus_protein_blasting_table, V1 == protein)
  if (min(subset_protein$V11) < 1e-4){
    homology = "True"
  } else {
    homology = "False"
  }
  Bos_taurus_protein_homology <- rbind(Bos_taurus_protein_homology, c(protein, homology))
}
colnames(Bos_taurus_protein_homology) <- c("GenBank.Accession", "Bos.taurus.Homology")

##Homo sapiens
Homo_sapiens_protein_blasting_table <- read.table("UK1_Homo_sapiens_tabular_070323.txt")
Homo_sapiens_protein_homology = data.frame()
for (protein in unique(Homo_sapiens_protein_blasting_table$V1)) {
  subset_protein = filter(Homo_sapiens_protein_blasting_table, V1 == protein)
  if (min(subset_protein$V11) < 1e-4){
    homology = "True"
  } else {
    homology = "False"
  }
  Homo_sapiens_protein_homology <- rbind(Homo_sapiens_protein_homology, c(protein, homology))
}
colnames(Homo_sapiens_protein_homology) <- c("GenBank.Accession", "Homo.sapiens.Homology")

#Merge the results from each of the potential hosts by GenBank.Accession columns
#all=TRUE allows for all results to be kept and missing information to be filled by NAs
negative_protein_homology = merge(Homo_sapiens_protein_homology, Mus_musculus_protein_homology, by="GenBank.Accession", all= TRUE)
negative_protein_homology = merge(negative_protein_homology, Sus_scrofa_protein_homology, by="GenBank.Accession", all= TRUE)
negative_protein_homology = merge(negative_protein_homology, Gallus_gallus_protein_homology, by="GenBank.Accession", all= TRUE)
negative_protein_homology = merge(negative_protein_homology, Meleagris_gallopavo_protein_homology, by="GenBank.Accession", all= TRUE)
negative_protein_homology = merge(negative_protein_homology, Bos_taurus_protein_homology, by="GenBank.Accession", all= TRUE)

#Change anything with NA to FALSE since that means the protein query did not have an evalue <10 (BLASTP Default cutoff)
negative_protein_homology[is.na(negative_protein_homology)] = "False"

####Positive homology####

##Make a table summarizing homology to Infantis
Infantis_protein_blasting_table <- read.table("UK1_Infantis_tabular_070323.txt")            #Upload your tabular results 
Infantis_protein_homology = data.frame()                                               #Create a blank dataframe to store summarized results
for (protein in unique(Infantis_protein_blasting_table$V1)) {                                   #For every unique protein query (UK1 Accession)...
  subset_protein = filter(Infantis_protein_blasting_table, V1 == protein)                       #Create a subset of the original dataframe with just that protein query
  if (max(subset_protein$V3) > 98){                                                #If the maximum protein identity is <98 then...
    homology = "True"                                                              #State "True" aka protein is homologous to serovar proteins
  } else {
    homology = "False"                                                             #Otherwise state "False" aka protein is nonhomologous to host proteins
  }
  Infantis_protein_homology <- rbind(Infantis_protein_homology, c(protein, homology))      #Save protein query and summary as a line in summarized dataframe
}
colnames(Infantis_protein_homology) <- c("GenBank.Accession", "Infantis.Homology")     #Change column names to something readable


##Make a table summarizing homology to Kentucky
Kentucky_protein_blasting_table <- read.table("UK1_Kentucky_tabular_070323.txt")
Kentucky_protein_homology = data.frame()
for (protein in unique(Kentucky_protein_blasting_table$V1)) {
  subset_protein = filter(Kentucky_protein_blasting_table, V1 == protein)
  if (max(subset_protein$V3) > 98){
    homology = "True"
  } else {
    homology = "False"
  }
  Kentucky_protein_homology <- rbind(Kentucky_protein_homology, c(protein, homology))
}
colnames(Kentucky_protein_homology) <- c("GenBank.Accession", "Kentucky.Homology")


##Make a table summarizing homology to Hadar
Hadar_protein_blasting_table <- read.table("UK1_Hadar_tabular_070323.txt")
Hadar_protein_homology = data.frame()
for (protein in unique(Hadar_protein_blasting_table$V1)) {
  subset_protein = filter(Hadar_protein_blasting_table, V1 == protein)
  if (max(subset_protein$V3) > 98){
    homology = "True"
  } else {
    homology = "False"
  }
  Hadar_protein_homology <- rbind(Hadar_protein_homology, c(protein, homology))
}
colnames(Hadar_protein_homology) <- c("GenBank.Accession", "Hadar.Homology")

##Make a table summarizing homology to Enteritidis
Enteritidis_protein_blasting_table <- read.table("UK1_Enteritidis_tabular_070323.txt")
Enteritidis_protein_homology = data.frame()
for (protein in unique(Enteritidis_protein_blasting_table$V1)) {
  subset_protein = filter(Enteritidis_protein_blasting_table, V1 == protein)
  if (max(subset_protein$V3) > 98){
    homology = "True"
  } else {
    homology = "False"
  }
  Enteritidis_protein_homology <- rbind(Enteritidis_protein_homology, c(protein, homology))
}
colnames(Enteritidis_protein_homology) <- c("GenBank.Accession", "Enteritidis.Homology")

##Make a table summarizing homology to Uganda
Uganda_protein_blasting_table <- read.table("UK1_Uganda_tabular_070323.txt")
Uganda_protein_homology = data.frame()
for (protein in unique(Uganda_protein_blasting_table$V1)) {
  subset_protein = filter(Uganda_protein_blasting_table, V1 == protein)
  if (max(subset_protein$V3) > 98){
    homology = "True"
  } else {
    homology = "False"
  }
  Uganda_protein_homology <- rbind(Uganda_protein_homology, c(protein, homology))
}
colnames(Uganda_protein_homology) <- c("GenBank.Accession", "Uganda.Homology")

#Merge the results from each of the potential hosts by GenBank.Accession columns
positive_protein_homology = merge(Enteritidis_protein_homology, Hadar_protein_homology, by="GenBank.Accession", all= TRUE)
positive_protein_homology = merge(positive_protein_homology, Infantis_protein_homology, by="GenBank.Accession", all= TRUE)
positive_protein_homology = merge(positive_protein_homology, Kentucky_protein_homology, by="GenBank.Accession", all= TRUE)
positive_protein_homology = merge(positive_protein_homology, Uganda_protein_homology, by="GenBank.Accession", all= TRUE)

#Change anything with NA to FALSE since that means the protein query did not have an evalue <10 (BLASTP Default cutoff)
positive_protein_homology[is.na(positive_protein_homology)] = "False"

##Merge positive and negative homology
full_protein_homology <- merge(negative_protein_homology, positive_protein_homology, by="GenBank.Accession", all = TRUE)

#Change anything with NA to FALSE since that means the protein query did not have an evalue <10 (BLASTP Default cutoff)
#This is especially important for the negative homology results where not every protein will have results,
#Whereas the positive homology results will have a much greater chance of having every protein represented
full_protein_homology[is.na(full_protein_homology)] = "False"

#Add in a row summarizing the results
full_protein_homology$Negative.Homology <- with(full_protein_homology, ifelse(Homo.sapiens.Homology == "True" | 
                                                                                Mus.musculus.Homology == "True" | 
                                                                                Sus.scrofa.Homology == "True" | 
                                                                                Gallus.gallus.Homology == "True" | 
                                                                                Meleagris.gallopavo.Homology =="True" | 
                                                                                Bos.taurus.Homology == "True",
                                                                              "True", "False"))

full_protein_homology$Positive.Homology <- with(full_protein_homology, ifelse(Enteritidis.Homology == "False" | 
                                                                                Hadar.Homology == "False" | 
                                                                                Infantis.Homology == "False" | 
                                                                                Kentucky.Homology == "False" | 
                                                                                Uganda.Homology =="False",
                                                                              "False", "True"))



#Copy results and add stepwise to Excel Sheet (i.e. do negative homology, then positive homology)
#Then see Reverse_Vaccinology_UK1_Scripts.txt for information regarding metadata columns to add
clipr::write_clip(full_protein_homology)

#...

#Importing information from Excel/Notepad++...
###Non-Similar to Hosts Proteins = 110###
###Similar to Other Salmonella Serovar Proteins = 101###

#...

####Adding Additional Annotations####

#Upload the Post-Filtering Proteins to create new metadata columns and add in additional annotations
Post_Filtering_UK1_Proteins <- read.delim("Post_Filtering_UK1_Proteins.txt")

##Get RefSeq Annotations and GO Terms for future reference (Construct Creation)

##Create a dataframe of exact sequence matches between GenBank and RefSeq

#Upload the Protein FASTA (.faa) files from the NCBI_Resources
Genbank <- readAAStringSet("../NCBI_Resources/GenBank_GCA_000213635.1_ASM21363v1/GCA_000213635.1_ASM21363v1_protein.faa")
RefSeq <- readAAStringSet("../NCBI_Resources/RefSeq_GCF_000213635.1_ASM21363v1/GCF_000213635.1_ASM21363v1_protein.faa")

#Convert to a dataframe and rename columns something useful
Genbank <- data.frame(names(Genbank), paste(Genbank))
colnames(Genbank) <- c("GenBank.Description", "Protein.Sequence")
RefSeq <- data.frame(names(RefSeq), paste(RefSeq))
colnames(RefSeq) <- c("RefSeq.Description", "Protein.Sequence")

#Merge the dataframes by sequence to find same annotations for exact sequences
Genbank_RefSeq <- merge(Genbank, RefSeq, by="Protein.Sequence")

#Separate resulting columns into Accession and Description 
Genbank_RefSeq <- separate(data = Genbank_RefSeq, col = GenBank.Description, into = c("GenBank.Accession", "GenBank.Protein.Description"), sep = " ", extra = "merge")
Genbank_RefSeq <- separate(data = Genbank_RefSeq, col = RefSeq.Description, into = c("RefSeq.Accession", "RefSeq.Protein.Description"), sep = " ", extra = "merge")

##Create a dataframe with GO terms from RefSeq

#Upload the RefSeq .gff file
RefSeq_gff <- read.gff("../NCBI_Resources/RefSeq_GCF_000213635.1_ASM21363v1/GCF_000213635.1_ASM21363v1_genomic.gff")

#Keep only the columns with Protein Homology information
RefSeq_gff_PH <- filter(RefSeq_gff, source == "Protein Homology")

##Make a dataframe of the GO terms

#Make a dataframe with the information in the attributes splite into two columns by "Ontology_term"
RefSeq_Ontology_terms <- as.data.frame(str_split_fixed(RefSeq_gff_PH$attributes, "Ontology_term=", 2))

#Keep the first element in the first column split by ";"
RefSeq_Ontology_terms$V1 <- sapply(strsplit(RefSeq_Ontology_terms$V1,";"), `[`, 1)

#Keep the 2nd element of the first column split by "-"
RefSeq_Ontology_terms$V1 <- sapply(strsplit(RefSeq_Ontology_terms$V1,"-"), `[`, 2)

#Keep the 1st element of the 2nd column split by ";"
RefSeq_Ontology_terms$V2 <- sapply(strsplit(RefSeq_Ontology_terms$V2,";"), `[`, 1)

#Rename columns something useful
colnames(RefSeq_Ontology_terms) <- c("RefSeq.Accession", "Ontology.Terms")

##Make a dataframe of GO functions

#Make a dataframe with the information in the attributes splite into two columns by "go_function"
RefSeq_GO_Functions <- as.data.frame(str_split_fixed(RefSeq_gff_PH$attributes, "go_function=", 2))

#Keep the 1st element of the first column split by ";"
RefSeq_GO_Functions$V1 <- sapply(strsplit(RefSeq_GO_Functions$V1,";"), `[`, 1)

#Keep the 2nd element of the first column split by "-"
RefSeq_GO_Functions$V1 <- sapply(strsplit(RefSeq_GO_Functions$V1,"-"), `[`, 2)

#Keep the 1st element of the 2nd column split by ";"
RefSeq_GO_Functions$V2 <- sapply(strsplit(RefSeq_GO_Functions$V2,";"), `[`, 1)

#Change the column names to something useful
colnames(RefSeq_GO_Functions) <- c("RefSeq.Accession", "GO.Functions")

#Merge the two GO associated dataframes
RefSeq_GO_Annotations <- merge(RefSeq_Ontology_terms, RefSeq_GO_Functions, by="RefSeq.Accession")

#Merge the combo GenBank/RefSeq dataframe wtih the GO dataframe
GenBank_RefSeq_GO_Annotations <- merge(Genbank_RefSeq, RefSeq_GO_Annotations, by="RefSeq.Accession")

#Change any NAs to " -"
GenBank_RefSeq_GO_Annotations[is.na(GenBank_RefSeq_GO_Annotations)] = " -"

#Merge with Post_Filtering_UK1_Proteins_Annotated while removing repeat columns
Post_Filtering_UK1_Proteins_Annotated <- merge(Post_Filtering_UK1_Proteins_Annotated, select(GenBank_RefSeq_GO_Annotations, c("GenBank.Accession", "RefSeq.Accession", "RefSeq.Protein.Description", "Ontology.Terms", "GO.Functions")), by="GenBank.Accession", all.x=TRUE)

#Change Order to put these new columns at the end
#Get original table order
Post_Annotation_Column_Order <- dput(colnames(Post_Filtering_UK1_Proteins))

#Make a list of new column names
Annotation_Column_Names <- c("Full.Category", "Gene.Name.Essentiality", "Localization.Essentiality", "Gene.Name.Localization", "RefSeq.Accession", "RefSeq.Protein.Description", "Ontology.Terms", "GO.Functions")

#Append the above columns names to the end of this list
Post_Annotation_Column_Order <- append(Post_Annotation_Column_Order, Annotation_Column_Names)

#Reorder columns
Post_Filtering_UK1_Proteins_Annotated <- Post_Filtering_UK1_Proteins_Annotated[, Post_Annotation_Column_Order] # leave the row index blank to keep all rows

#Copy down into Excel
clipr::write_clip(Post_Filtering_UK1_Proteins_Annotated)

####PHASE II EPITOPE FILTERING####

#Change working directory
setwd("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Projects_Reverse_Vaccinology/UK1/cleaned_up/Phase_II_Epitope_Discovery")

#...

#Importing information from Notepad++/Excel...
###855 Non-Unique MHCII Epitopes###

###885 Non-Unique MHCI Epitopes###
###396 Non-Unique Immunogenic MHCI Epitopes###

#...


####Antigencity####

#Create fasta versions of epitopes

#MHCI
Pre_VaxiJen_UK1_MHCI_Epitope_Peptides <- read.delim("Pre_VaxiJen_UK1_MHCI_Epitope_Peptides.txt ", header=FALSE)
Pre_VaxiJen_UK1_MHCI_Epitope_Unique_Ids <- read.delim("Pre_VaxiJen_UK1_MHCI_Epitope_Unique_Ids.txt", header=FALSE)
write.fasta(as.list(dput(Pre_VaxiJen_UK1_MHCI_Epitope_Peptides$V1)), as.list(dput(Pre_VaxiJen_UK1_MHCI_Epitope_Unique_Ids$V1)), "Pre_VaxiJen_UK1_MHCI_Epitopes.fasta")

#MHCII
Pre_VaxiJen_UK1_MHCII_Epitope_Peptides <- read.delim("Pre_VaxiJen_UK1_MHCII_Epitope_Peptides.txt ", header=FALSE)
Pre_VaxiJen_UK1_MHCII_Epitope_Unique_Ids <- read.delim("Pre_VaxiJen_UK1_MHCII_Epitope_Unique_Ids.txt", header=FALSE)
write.fasta(as.list(dput(Pre_VaxiJen_UK1_MHCII_Epitope_Peptides$V1)), as.list(dput(Pre_VaxiJen_UK1_MHCII_Epitope_Unique_Ids$V1)), "Pre_VaxiJen_UK1_MHCII_Epitopes.fasta")

#Head back over to Reverse_Vaccinology_UK1_Scripts.txt

#...

#Importing information from Notepad++/Excel...
###210 Non-Unique MHCI Antigenic Epitopes###
###210 Non-Unique MHCI Non-Toxic Epitopes###
###132 Non-Unique MHCI Hydrophilic Epitopes###

###446 Non-Unique MHCII Antigenic Epitopes###
###440 Non-Unique MHCII Non-Toxic Epitopes###
###394 Non-Unique MHCII Hydrophilic Epitopes###

#...

####Positive Blasting####

#...

#Since this blasting will be a combination of MHCI and MHCII epitopes and we are going to need this data anyways for BepiPred,
#save the Pre Positive Homology Results as tab-deliminated txts and upload into R
Pre_Positive_Homology_UK1_MHCI_Results <- read.delim("Pre_Positive_Homology_UK1_MHCI_Results.txt", header = TRUE)
Pre_Positive_Homology_UK1_MHCII_Results <- read.delim("Pre_Positive_Homology_UK1_MHCII_Results.txt", header = TRUE)

#Create fastas from these results to run use while blasting

#MHCI
write.fasta(as.list(dput(Pre_Positive_Homology_UK1_MHCI_Results$Peptide)), as.list(dput(Pre_Positive_Homology_UK1_MHCI_Results$Unique_Id)), "Pre_Positive_Homology_UK1_MHCI_Epitopes.fasta")

#MHCII
write.fasta(as.list(dput(Pre_Positive_Homology_UK1_MHCII_Results$Peptide)), as.list(dput(Pre_Positive_Homology_UK1_MHCII_Results$Unique_Id)), "Pre_Positive_Homology_UK1_MHCII_Epitopes.fasta")

#Head back over to Reverse_Vaccinology_UK1_Scripts.txt

#...

#Change working directory
setwd("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Projects_Reverse_Vaccinology/UK1/cleaned_up/Phase_II_Epitope_Discovery/epitope_proteome_blasting")

#Upload <serovar>_proteomes_epitope_homology_summary.txt files from SCINet 
Enterititis_proteomes_epitope_homology_summary <- read.delim("Enterititis/Enterititis_proteomes_epitope_homology_summary.txt", sep = " ")
Hadar_proteomes_epitope_homology_summary <- read.delim("Hadar/Hadar_proteomes_epitope_homology_summary.txt", sep = " ")
Infantis_proteomes_epitope_homology_summary <- read.delim("Infantis/Infantis_proteomes_epitope_homology_summary.txt", sep = " ")
Kentucky_proteomes_epitope_homology_summary <- read.delim("Kentucky/Kentucky_proteomes_epitope_homology_summary.txt", sep = " ")
Typhimurium_proteomes_epitope_homology_summary <- read.delim("Typhimurium/Typhimurium_proteomes_epitope_homology_summary.txt", sep = " ")
Uganda_proteomes_epitope_homology_summary <- read.delim("Uganda/Uganda_proteomes_epitope_homology_summary.txt", sep = " ")



#Make histograms to view the distribution of 100% homology across proteomes across epitopes (All %s)
Enterititis_epitopes_vs_proteomes_histo <- ggplot(Enteritidis_proteomes_epitope_homology_summary, aes(x=Epitope_vs_Enteritidis_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Enterititis_epitopes_vs_proteomes_histo

Hadar_epitopes_vs_proteomes_histo <- ggplot(Hadar_proteomes_epitope_homology_summary, aes(x=Epitope_vs_Hadar_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Hadar_epitopes_vs_proteomes_histo

Infantis_epitopes_vs_proteomes_histo <- ggplot(Infantis_proteomes_epitope_homology_summary, aes(x=Epitope_vs_Infantis_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Infantis_epitopes_vs_proteomes_histo

Kentucky_epitopes_vs_proteomes_histo <- ggplot(Kentucky_proteomes_epitope_homology_summary, aes(x=Epitope_vs_Kentucky_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Kentucky_epitopes_vs_proteomes_histo

Typhymiurium_epitopes_vs_proteomes_histo <- ggplot(Typhimurium_proteomes_epitope_homology_summary, aes(x=Epitope_vs_Typhimurium_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Typhymiurium_epitopes_vs_proteomes_histo

Uganda_epitopes_vs_proteomes_histo <- ggplot(Uganda_proteomes_epitope_homology_summary, aes(x=Epitope_vs_Uganda_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Uganda_epitopes_vs_proteomes_histo

#Make histograms to view the distribution of 100% homology across proteomes across epitopes (>90%)
Enteritidis_epitopes_vs_proteomes_histo_95 <- ggplot(filter(Enteritidis_proteomes_epitope_homology_summary, Epitope_vs_Enteritidis_proteomes_Homology_Percentage >95), aes(x=Epitope_vs_Enteritidis_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Enteritidis_epitopes_vs_proteomes_histo_95

Hadar_epitopes_vs_proteomes_histo_95 <- ggplot(filter(Hadar_proteomes_epitope_homology_summary, Epitope_vs_Hadar_proteomes_Homology_Percentage >95), aes(x=Epitope_vs_Hadar_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Hadar_epitopes_vs_proteomes_histo_95

Infantis_epitopes_vs_proteomes_histo_95 <- ggplot(filter(Infantis_proteomes_epitope_homology_summary, Epitope_vs_Infantis_proteomes_Homology_Percentage >95), aes(x=Epitope_vs_Infantis_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Infantis_epitopes_vs_proteomes_histo_95

Kentucky_epitopes_vs_proteomes_histo_95 <- ggplot(filter(Kentucky_proteomes_epitope_homology_summary, Epitope_vs_Kentucky_proteomes_Homology_Percentage >95), aes(x=Epitope_vs_Kentucky_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Kentucky_epitopes_vs_proteomes_histo_95

Typhimurium_epitopes_vs_proteomes_histo_95 <- ggplot(filter(Typhimurium_proteomes_epitope_homology_summary, Epitope_vs_Typhimurium_proteomes_Homology_Percentage >95), aes(x=Epitope_vs_Typhimurium_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Typhimurium_epitopes_vs_proteomes_histo_95

Uganda_epitopes_vs_proteomes_histo_95 <- ggplot(filter(Uganda_proteomes_epitope_homology_summary, Epitope_vs_Uganda_proteomes_Homology_Percentage >95), aes(x=Epitope_vs_Uganda_proteomes_Homology_Percentage)) +
  geom_histogram() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
Uganda_epitopes_vs_proteomes_histo_95

#Get the distributions of epitopes that did or did not pass the various proteome percentages filters for each serovar
Enteritidis_proteomes_epitope_homology_counts_summary <- as.data.frame(cbind(t(table(Enteritidis_proteomes_epitope_homology_summary$Epitope_vs_Enteritidis_proteomes_100_percent_Homology)), 
                                                                             t(table(Enteritidis_proteomes_epitope_homology_summary$Epitope_vs_Enteritidis_proteomes_99_percent_Homology)), 
                                                                             t(table(Enteritidis_proteomes_epitope_homology_summary$Epitope_vs_Enteritidis_proteomes_98_percent_Homology)), 
                                                                             t(table(Enteritidis_proteomes_epitope_homology_summary$Epitope_vs_Enteritidis_proteomes_97_percent_Homology)), 
                                                                             t(table(Enteritidis_proteomes_epitope_homology_summary$Epitope_vs_Enteritidis_proteomes_96_percent_Homology)), 
                                                                             t(table(Enteritidis_proteomes_epitope_homology_summary$Epitope_vs_Enteritidis_proteomes_95_percent_Homology)), 
                                                                             t(table(Enteritidis_proteomes_epitope_homology_summary$Epitope_vs_Enteritidis_proteomes_90_percent_Homology))))


Hadar_proteomes_epitope_homology_counts_summary <- as.data.frame(cbind(t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_100_percent_Homology)), 
                                                                       t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_99_percent_Homology)), 
                                                                       t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_98_percent_Homology)), 
                                                                       t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_97_percent_Homology)), 
                                                                       t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_96_percent_Homology)), 
                                                                       t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_95_percent_Homology)), 
                                                                       t(table(Hadar_proteomes_epitope_homology_summary$Epitope_vs_Hadar_proteomes_90_percent_Homology))))

Infantis_proteomes_epitope_homology_counts_summary <- as.data.frame(cbind(t(table(Infantis_proteomes_epitope_homology_summary$Epitope_vs_Infantis_proteomes_100_percent_Homology)), 
                                                                          t(table(Infantis_proteomes_epitope_homology_summary$Epitope_vs_Infantis_proteomes_99_percent_Homology)), 
                                                                          t(table(Infantis_proteomes_epitope_homology_summary$Epitope_vs_Infantis_proteomes_98_percent_Homology)), 
                                                                          t(table(Infantis_proteomes_epitope_homology_summary$Epitope_vs_Infantis_proteomes_97_percent_Homology)), 
                                                                          t(table(Infantis_proteomes_epitope_homology_summary$Epitope_vs_Infantis_proteomes_96_percent_Homology)), 
                                                                          t(table(Infantis_proteomes_epitope_homology_summary$Epitope_vs_Infantis_proteomes_95_percent_Homology)), 
                                                                          t(table(Infantis_proteomes_epitope_homology_summary$Epitope_vs_Infantis_proteomes_90_percent_Homology))))

Kentucky_proteomes_epitope_homology_counts_summary <- as.data.frame(cbind(t(table(Kentucky_proteomes_epitope_homology_summary$Epitope_vs_Kentucky_proteomes_100_percent_Homology)), 
                                                                          t(table(Kentucky_proteomes_epitope_homology_summary$Epitope_vs_Kentucky_proteomes_99_percent_Homology)), 
                                                                          t(table(Kentucky_proteomes_epitope_homology_summary$Epitope_vs_Kentucky_proteomes_98_percent_Homology)), 
                                                                          t(table(Kentucky_proteomes_epitope_homology_summary$Epitope_vs_Kentucky_proteomes_97_percent_Homology)), 
                                                                          t(table(Kentucky_proteomes_epitope_homology_summary$Epitope_vs_Kentucky_proteomes_96_percent_Homology)), 
                                                                          t(table(Kentucky_proteomes_epitope_homology_summary$Epitope_vs_Kentucky_proteomes_95_percent_Homology)), 
                                                                          t(table(Kentucky_proteomes_epitope_homology_summary$Epitope_vs_Kentucky_proteomes_90_percent_Homology))))

Typhimurium_proteomes_epitope_homology_counts_summary <- as.data.frame(cbind(t(table(Typhimurium_proteomes_epitope_homology_summary$Epitope_vs_Typhimurium_proteomes_100_percent_Homology)), 
                                                                             t(table(Typhimurium_proteomes_epitope_homology_summary$Epitope_vs_Typhimurium_proteomes_99_percent_Homology)), 
                                                                             t(table(Typhimurium_proteomes_epitope_homology_summary$Epitope_vs_Typhimurium_proteomes_98_percent_Homology)), 
                                                                             t(table(Typhimurium_proteomes_epitope_homology_summary$Epitope_vs_Typhimurium_proteomes_97_percent_Homology)), 
                                                                             t(table(Typhimurium_proteomes_epitope_homology_summary$Epitope_vs_Typhimurium_proteomes_96_percent_Homology)), 
                                                                             t(table(Typhimurium_proteomes_epitope_homology_summary$Epitope_vs_Typhimurium_proteomes_95_percent_Homology)), 
                                                                             t(table(Typhimurium_proteomes_epitope_homology_summary$Epitope_vs_Typhimurium_proteomes_90_percent_Homology))))

Uganda_proteomes_epitope_homology_counts_summary <- as.data.frame(cbind(t(table(Uganda_proteomes_epitope_homology_summary$Epitope_vs_Uganda_proteomes_100_percent_Homology)), 
                                                                        t(table(Uganda_proteomes_epitope_homology_summary$Epitope_vs_Uganda_proteomes_99_percent_Homology)), 
                                                                        t(table(Uganda_proteomes_epitope_homology_summary$Epitope_vs_Uganda_proteomes_98_percent_Homology)), 
                                                                        t(table(Uganda_proteomes_epitope_homology_summary$Epitope_vs_Uganda_proteomes_97_percent_Homology)), 
                                                                        t(table(Uganda_proteomes_epitope_homology_summary$Epitope_vs_Uganda_proteomes_96_percent_Homology)), 
                                                                        t(table(Uganda_proteomes_epitope_homology_summary$Epitope_vs_Uganda_proteomes_95_percent_Homology)), 
                                                                        t(table(Uganda_proteomes_epitope_homology_summary$Epitope_vs_Uganda_proteomes_90_percent_Homology))))



#Compare and contrast the different levels of proteome homology across all serovars at once
test_positive_epitope_homology <- merge(Kentucky_proteomes_epitope_homology_summary, Hadar_proteomes_epitope_homology_summary)
test_positive_epitope_homology <- merge(test_positive_epitope_homology, Uganda_proteomes_epitope_homology_summary)
test_positive_epitope_homology <- merge(test_positive_epitope_homology, Infantis_proteomes_epitope_homology_summary)
test_positive_epitope_homology <- merge(test_positive_epitope_homology, Typhimurium_proteomes_epitope_homology_summary)
test_positive_epitope_homology <- merge(test_positive_epitope_homology, Enteritidis_proteomes_epitope_homology_summary)

#Explore the distribution of TRUE and FALSE epitopes

#100%
test_positive_epitope_homology$Positive_Epitope_Homology_100 <- with(test_positive_epitope_homology, ifelse(Epitope_vs_Kentucky_proteomes_100_percent_Homology == "FALSE" | 
                                                                                                              Epitope_vs_Hadar_proteomes_100_percent_Homology == "FALSE" |
                                                                                                              Epitope_vs_Uganda_proteomes_100_percent_Homology == "FALSE" |
                                                                                                              Epitope_vs_Infantis_proteomes_100_percent_Homology == "FALSE" |
                                                                                                              Epitope_vs_Infantis_proteomes_100_percent_Homology == "FALSE" |
                                                                                                              Epitope_vs_Typhimurium_proteomes_100_percent_Homology == "FALSE" |
                                                                                                              Epitope_vs_Enteritidis_proteomes_100_percent_Homology == "FALSE",
                                                                                                            "FALSE", "TRUE"))
table(test_positive_epitope_homology$Positive_Epitope_Homology_100)
# FALSE  TRUE 
# 526     0

#99%
test_positive_epitope_homology$Positive_Epitope_Homology_99 <- with(test_positive_epitope_homology, ifelse(Epitope_vs_Kentucky_proteomes_99_percent_Homology == "FALSE" | 
                                                                                                             Epitope_vs_Hadar_proteomes_99_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Uganda_proteomes_99_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Infantis_proteomes_99_percent_Homology == "FALSE"|
                                                                                                             Epitope_vs_Typhimurium_proteomes_99_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Enteritidis_proteomes_99_percent_Homology == "FALSE",
                                                                                                           "FALSE", "TRUE"))
table(test_positive_epitope_homology$Positive_Epitope_Homology_99)
# FALSE  TRUE 
# 241   285

#98%
test_positive_epitope_homology$Positive_Epitope_Homology_98 <- with(test_positive_epitope_homology, ifelse(Epitope_vs_Kentucky_proteomes_98_percent_Homology == "FALSE" | 
                                                                                                             Epitope_vs_Hadar_proteomes_98_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Uganda_proteomes_98_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Infantis_proteomes_98_percent_Homology == "FALSE"|
                                                                                                             Epitope_vs_Typhimurium_proteomes_98_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Enteritidis_proteomes_98_percent_Homology == "FALSE",
                                                                                                           "FALSE", "TRUE"))
table(test_positive_epitope_homology$Positive_Epitope_Homology_98)
# FALSE  TRUE 
# 196   330

#97%
test_positive_epitope_homology$Positive_Epitope_Homology_97 <- with(test_positive_epitope_homology, ifelse(Epitope_vs_Kentucky_proteomes_97_percent_Homology == "FALSE" | 
                                                                                                             Epitope_vs_Hadar_proteomes_97_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Uganda_proteomes_97_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Infantis_proteomes_97_percent_Homology == "FALSE"|
                                                                                                             Epitope_vs_Typhimurium_proteomes_97_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Enteritidis_proteomes_97_percent_Homology == "FALSE",
                                                                                                           "FALSE", "TRUE"))
table(test_positive_epitope_homology$Positive_Epitope_Homology_97)
# FALSE  TRUE 
# 163   363

#96%
test_positive_epitope_homology$Positive_Epitope_Homology_96 <- with(test_positive_epitope_homology, ifelse(Epitope_vs_Kentucky_proteomes_96_percent_Homology == "FALSE" | 
                                                                                                             Epitope_vs_Hadar_proteomes_96_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Uganda_proteomes_96_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Infantis_proteomes_96_percent_Homology == "FALSE"|
                                                                                                             Epitope_vs_Typhimurium_proteomes_96_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Enteritidis_proteomes_96_percent_Homology == "FALSE",
                                                                                                           "FALSE", "TRUE"))
table(test_positive_epitope_homology$Positive_Epitope_Homology_96)
# FALSE  TRUE 
# 136   390

#95%
test_positive_epitope_homology$Positive_Epitope_Homology_95 <- with(test_positive_epitope_homology, ifelse(Epitope_vs_Kentucky_proteomes_95_percent_Homology == "FALSE" | 
                                                                                                             Epitope_vs_Hadar_proteomes_95_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Uganda_proteomes_95_percent_Homology == "FALSE"|
                                                                                                             Epitope_vs_Infantis_proteomes_95_percent_Homology == "FALSE"|
                                                                                                             Epitope_vs_Typhimurium_proteomes_95_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Enteritidis_proteomes_95_percent_Homology == "FALSE",
                                                                                                           "FALSE", "TRUE"))
table(test_positive_epitope_homology$Positive_Epitope_Homology_95)
# FALSE  TRUE 
# 135   391

#90%
test_positive_epitope_homology$Positive_Epitope_Homology_90 <- with(test_positive_epitope_homology, ifelse(Epitope_vs_Kentucky_proteomes_90_percent_Homology == "FALSE" | 
                                                                                                             Epitope_vs_Hadar_proteomes_90_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Uganda_proteomes_90_percent_Homology == "FALSE"|
                                                                                                             Epitope_vs_Infantis_proteomes_90_percent_Homology == "FALSE"|
                                                                                                             Epitope_vs_Typhimurium_proteomes_90_percent_Homology == "FALSE" |
                                                                                                             Epitope_vs_Enteritidis_proteomes_90_percent_Homology == "FALSE",
                                                                                                           "FALSE", "TRUE"))
table(test_positive_epitope_homology$Positive_Epitope_Homology_90)
# FALSE  TRUE 
# 132   394

#Based upon the results from all 6 serovars, the best combination of surviving epitopes and stringent filters is >=99% of proteomes with 100% identity to an epitope

#Keep the relevant columns
positive_epitope_homology <- select(test_positive_epitope_homology, c("Unique_Id", "Identity", 
                                                                      "Epitope_vs_Kentucky_proteomes_Homology_Percentage", "Epitope_vs_Kentucky_proteomes_99_percent_Homology", 
                                                                      "Epitope_vs_Hadar_proteomes_Homology_Percentage", "Epitope_vs_Hadar_proteomes_99_percent_Homology", 
                                                                      "Epitope_vs_Uganda_proteomes_Homology_Percentage", "Epitope_vs_Uganda_proteomes_99_percent_Homology", 
                                                                      "Epitope_vs_Infantis_proteomes_Homology_Percentage", "Epitope_vs_Infantis_proteomes_99_percent_Homology", 
                                                                      "Epitope_vs_Typhimurium_proteomes_Homology_Percentage", "Epitope_vs_Typhimurium_proteomes_99_percent_Homology", 
                                                                      "Epitope_vs_Enteritidis_proteomes_Homology_Percentage", "Epitope_vs_Enteritidis_proteomes_99_percent_Homology"))

#Change anything with NA to FALSE since that means the protein query did not have an evalue <10000 (BLASTP Default cutoff)
positive_epitope_homology[is.na(positive_epitope_homology)] = "FALSE"

#Add in the final column summarizing the results
positive_epitope_homology$Positive_Epitope_Homology <- with(positive_epitope_homology, ifelse(Epitope_vs_Kentucky_proteomes_99_percent_Homology == "FALSE" | 
                                                                                                Epitope_vs_Hadar_proteomes_99_percent_Homology == "FALSE" |
                                                                                                Epitope_vs_Uganda_proteomes_99_percent_Homology == "FALSE" |
                                                                                                Epitope_vs_Infantis_proteomes_99_percent_Homology == "FALSE" |
                                                                                                Epitope_vs_Typhimurium_proteomes_99_percent_Homology == "FALSE" |
                                                                                                Epitope_vs_Enteritidis_proteomes_99_percent_Homology == "FALSE",
                                                                                              "FALSE", "TRUE"))

##MHCI
#Merge the positive homology results with the pre positive homology results
Post_Positive_Homology_UK1_MHCI_Results <- merge(Pre_Positive_Homology_UK1_MHCI_Results, positive_epitope_homology) #132 rows

#Keep epitopes with 100% percent identity to 99% of all serovars' proteomes
Post_Filtering_UK1_MHCI_Results <- filter(Post_Positive_Homology_UK1_MHCI_Results, Positive_Epitope_Homology =="TRUE") #69 rows

#Determine number of unique epitopes lost
length(unique(Post_Positive_Homology_UK1_MHCI_Results$Peptide)) #76
length(unique(Post_Filtering_UK1_MHCI_Results$Peptide)) #42

#Determine number of unique proteins lost
length(unique(Post_Positive_Homology_UK1_MHCI_Results$Identity)) #51
length(unique(Post_Filtering_UK1_MHCI_Results$Identity)) #34

##MHCII
#Merge the positive homology results with the pre positive homology results
Post_Positive_Homology_UK1_MHCII_Results <- merge(Pre_Positive_Homology_UK1_MHCII_Results, positive_epitope_homology) #394 rows

#Keep epitopes with 100% percent identity to other serovars
Post_Filtering_UK1_MHCII_Results <- filter(Post_Positive_Homology_UK1_MHCII_Results, Positive_Epitope_Homology =="TRUE") #256 rows

#Determine number of unique epitopes lost
length(unique(Post_Positive_Homology_UK1_MHCII_Results$Peptide)) #231
length(unique(Post_Filtering_UK1_MHCII_Results$Peptide)) #126

#Determine number of unique proteins lost
length(unique(Post_Positive_Homology_UK1_MHCII_Results$Identity)) #57
length(unique(Post_Filtering_UK1_MHCII_Results$Identity)) #41


#Transfer information of four above dataframes back down into Excel
clipr::write_clip(Post_Positive_Homology_UK1_MHCI_Results)
clipr::write_clip(Post_Positive_Homology_UK1_MHCII_Results)

#Merge with annotation information for analysis of proteins that do or do not make it
test_positive_epitope_homology_annotated <- merge(test_positive_epitope_homology, Post_Filtering_UK1_Proteins_Annotated_Epitope_Ver) #394 rows

#Make a dataframe of the epitopes that were false across all serovars
total_false_positive_epitope_homology_annotated <- filter(test_positive_epitope_homology_annotated, Epitope_vs_Kentucky_proteomes_90_percent_Homology == "FALSE" & 
                                                            Epitope_vs_Hadar_proteomes_90_percent_Homology == "FALSE" &
                                                            Epitope_vs_Uganda_proteomes_90_percent_Homology == "FALSE" &
                                                            Epitope_vs_Infantis_proteomes_90_percent_Homology == "FALSE" &
                                                            Epitope_vs_Typhimurium_proteomes_90_percent_Homology == "FALSE" &
                                                            Epitope_vs_Enteritidis_proteomes_90_percent_Homology == "FALSE")
#9

#Head back over to Reverse_Vaccinology_UK1_Scripts.txt for Linear B-Cell (LBL) epitope identification

#...

####BepiPred Linear B-Cell Epitope Identification####

##Proteins
UK1_bepipred_raw_results <- read.table("Bepipred_Raw_Results.txt", header = TRUE)             #Read in raw BepiPred results
UK1_bepipred_sequence_results = data.frame()                                                  #Create a blank dataframe to hold summarized results
for (protein in unique(UK1_bepipred_raw_results$Accession)) {                                 #For each unique protein in the BepiPred results...
  subset_protein = filter(UK1_bepipred_raw_results, Accession == protein)                     #Subset the dataframe to focus on that protein
  sequence <- ""                                                                              #Create a blank vector to store sequence results in
  for(row in 1:nrow(subset_protein)) {                                                        #For every row in this subset dataframe...
    residue <- subset_protein[row, "Residue"]                                                 #Save the residue
    score <- subset_protein[row, "Rolling.Mean.Score"]                                        #Save the score
    if(score>0.1512) {                                                                        #If the score is greater than the default threshold...
      residue <- toupper(residue)                                                             #Convert the residue to upper case
    } else {
      residue <- tolower(residue)                                                             #Otherwise, convert residue to lower case
    }
    sequence <- paste(sequence, residue, sep = "")                                            #Then paste the sequence together with no deliminator
  }
  UK1_bepipred_sequence_results <- rbind(UK1_bepipred_sequence_results, c(protein, sequence)) #Save the unique protein the the Bepipred sequence in the blank dataframe
}
colnames(UK1_bepipred_sequence_results) <- c("Identity", "Bepipred_Sequence")                 #Change the column names to something readable

#Remove intermediates
remove(protein)
remove(sequence)
remove(row)
remove(score)
remove(residue)
remove(subset_protein)

##MHCI Epitopes
#Merge the results with the post filtering dataframe since you are gonna need info from that (Will merge by Identity aka GenBank Accession)
Post_Filtering_UK1_MHCI_BepiPred <- merge(Post_Filtering_UK1_MHCI_Results, UK1_bepipred_sequence_results)

Post_Filtering_UK1_MHCI_BepiPred_Percent <-data.frame()                                            #Create a blank dataframe to store results
for(row in 1:nrow(Post_Filtering_UK1_MHCI_BepiPred)){                                              #For every row in the merged dataframe from above...
  test_epitope <- Post_Filtering_UK1_MHCI_BepiPred[row, "Peptide"]                                 #Save the Peptide as the test_epitope
  test_ptn <- Post_Filtering_UK1_MHCI_BepiPred[row, "Bepipred_Sequence"]                           #Save the BepiPred Sequence as the test_ptn
  test_unique_id <- Post_Filtering_UK1_MHCI_BepiPred[row, "Unique_Id"]                             #Save the  Unique_Id as the test_unique_id
  detected_epitope <- str_extract(test_ptn, fixed(test_epitope, ignore_case=TRUE))                 #Extract out portion of the BepiPred Sequence that matches the Peptide while keeping the informative letter case from the BepiPred Sequence
  Bcell_epitope_percentage <- str_count(detected_epitope, "[A-Z]") / nchar(detected_epitope) *100  #Determine what percentage of the extracted peptide is upper case and thus a projected Linear B-Cell Residue
  Post_Filtering_UK1_MHCI_BepiPred_Percent <- rbind(Post_Filtering_UK1_MHCI_BepiPred_Percent, c(test_unique_id, test_epitope, test_ptn, detected_epitope, Bcell_epitope_percentage)) #Save everything in a dataframe
}

#Change column names to something readable
colnames(Post_Filtering_UK1_MHCI_BepiPred_Percent) <- c("Unique_Id", "Peptide", "Bepipred_Sequence", "detected_epitope", "Bcell_epitope_percentage") 

#Merge the resulting dataframe with the post filtering results with the BepiPred Sequence while removing duplicate columns
Post_Filtering_UK1_MHCI_Full <- merge(Post_Filtering_UK1_MHCI_BepiPred, select(Post_Filtering_UK1_MHCI_BepiPred_Percent, c("Unique_Id", "detected_epitope", "Bcell_epitope_percentage")), by="Unique_Id")

#Save in Excel
clipr::write_clip(Post_Filtering_UK1_MHCI_Full)

##MCHII Epitopes
#Merge the results with the post filtering dataframe since you are gonna need info from that (Will merge by Identity aka GenBank Accession)
Post_Filtering_UK1_MHCII_BepiPred <- merge(Post_Filtering_UK1_MHCII_Results, UK1_bepipred_sequence_results)

Post_Filtering_UK1_MHCII_BepiPred_Percent <-data.frame()
for(row in 1:nrow(Post_Filtering_UK1_MHCII_BepiPred)){
  test_epitope <- Post_Filtering_UK1_MHCII_BepiPred[row, "Peptide"]
  test_ptn <- Post_Filtering_UK1_MHCII_BepiPred[row, "Bepipred_Sequence"]
  test_unique_id <- Post_Filtering_UK1_MHCII_BepiPred[row, "Unique_Id"]
  detected_epitope <- str_extract(test_ptn, fixed(test_epitope, ignore_case=TRUE))
  Bcell_epitope_percentage <- str_count(detected_epitope, "[A-Z]") / nchar(detected_epitope) *100
  Post_Filtering_UK1_MHCII_BepiPred_Percent <- rbind(Post_Filtering_UK1_MHCII_BepiPred_Percent, c(test_unique_id, test_epitope, test_ptn, detected_epitope, Bcell_epitope_percentage))
}

#Remove intermediates
remove(test_unique_id)
remove(test_epitope)
remove(test_ptn)
remove(detected_epitope)
remove(Bcell_epitope_percentage)

#Change column names to something readable
colnames(Post_Filtering_UK1_MHCII_BepiPred_Percent) <- c("Unique_Id", "Peptide", "Bepipred_Sequence", "detected_epitope", "Bcell_epitope_percentage")

#Merge the resulting dataframe with the post filtering results with the BepiPred Sequence while removing duplicate columns
Post_Filtering_UK1_MHCII_Full <- merge(Post_Filtering_UK1_MHCII_BepiPred, select(Post_Filtering_UK1_MHCII_BepiPred_Percent, c("Unique_Id", "detected_epitope", "Bcell_epitope_percentage")), by="Unique_Id")

#Save in Excel
clipr::write_clip(Post_Filtering_UK1_MHCII_Full)



####Find Unique Epitopes####

#Select for columns common between MHCI and MHCII tables
Common_MHCI_MHCII_Column_Names <- intersect(dput(colnames(Post_Filtering_UK1_MHCI_Full)), dput(colnames(Post_Filtering_UK1_MHCII_Full)))

##Create a MHCI table summarizing information per unique peptide
#Summarize the EL Score and Percent Rank for each unique peptide
Post_Filtering_UK1_MHCI_Unique_Peptide_Summary = data.frame()
for (protein in unique(Post_Filtering_UK1_MHCI_Full$Peptide)) {                   #For each unique protein...
  subset_protein = filter(Post_Filtering_UK1_MHCI_Full, Peptide == protein)       #Subset the main dataframe to just that protein
  avg_score_el <- mean(subset_protein$Score_EL)                                    #Get the average Score EL
  avg_percent_rank_el <- mean(subset_protein$X.Rank_EL)                            #Get the average Percent Rank EL
  Post_Filtering_UK1_MHCI_Unique_Peptide_Summary <- rbind(Post_Filtering_UK1_MHCI_Unique_Peptide_Summary, c(protein, avg_score_el, avg_percent_rank_el)) #Save the three elements as a dataframe
}

#Make the column names readable
colnames(Post_Filtering_UK1_MHCI_Unique_Peptide_Summary) <- c("Peptide", "Avg_Score_EL", "Avg_X.Rank_EL")

#Create a dataframe that summarizes all information by Peptide
#This will take unique items in each column and make them into a list
temp_df <- select(Post_Filtering_UK1_MHCI_Full, all_of(Common_MHCI_MHCII_Column_Names)) %>%
  group_by(Peptide) %>%
  # Combine all strings
  summarize_all(funs(toString(unique(.[!is.na(.)])))) %>%
  ungroup()

#Combine the two dataframes while removing repeats and Unique_Id which we will remake
Post_Filtering_UK1_MHCI_Unique_Peptide_Summary <- merge(select(temp_df, -c(Score_EL, X.Rank_EL, Unique_Id)), Post_Filtering_UK1_MHCI_Unique_Peptide_Summary)

#Remake the Unique_Id column with new combined MHC column
Post_Filtering_UK1_MHCI_Unique_Peptide_Summary$Unique_Id <- paste(Post_Filtering_UK1_MHCI_Unique_Peptide_Summary$Identity, "_MHCI_", Post_Filtering_UK1_MHCI_Unique_Peptide_Summary$MHC, "_Pos_", Post_Filtering_UK1_MHCI_Unique_Peptide_Summary$Pos, sep = "")

#Rename the MHC column MHC_Allele
names(Post_Filtering_UK1_MHCI_Unique_Peptide_Summary)[names(Post_Filtering_UK1_MHCI_Unique_Peptide_Summary) == 'MHC'] <- 'MHC_Allele'

#Make a new column called MHC_Type
Post_Filtering_UK1_MHCI_Unique_Peptide_Summary$MHC_Type <- "MHCI"

#Get the column names
dput(colnames(Post_Filtering_UK1_MHCI_Unique_Peptide_Summary))

#Put the column names in a readable order
Unique_Peptide_Summary_Column_Order <- c("Identity", "MHC_Type",  "MHC_Allele", "Pos", "Peptide", "Core", "Of",  
                                         "Avg_Score_EL", "Avg_X.Rank_EL", "BindLevel", "VaxiJen_Score", "Unique_Id", "ToxinPred_SVM_Score", "ToxinPred_Prediction", 
                                         "Hydrophobicity", "GRAVY_Hydropathicity", "Hydrophilicity", "Charge", 
                                         "Epitope_vs_Enteritidis_proteomes_Homology_Percentage", "Epitope_vs_Enteritidis_proteomes_99_percent_Homology", 
                                         "Epitope_vs_Hadar_proteomes_Homology_Percentage", "Epitope_vs_Hadar_proteomes_99_percent_Homology", 
                                         "Epitope_vs_Infantis_proteomes_Homology_Percentage", "Epitope_vs_Infantis_proteomes_99_percent_Homology", 
                                         "Epitope_vs_Kentucky_proteomes_Homology_Percentage", "Epitope_vs_Kentucky_proteomes_99_percent_Homology", 
                                         "Epitope_vs_Typhimurium_proteomes_Homology_Percentage", "Epitope_vs_Typhimurium_proteomes_99_percent_Homology", 
                                         "Epitope_vs_Uganda_proteomes_Homology_Percentage", "Epitope_vs_Uganda_proteomes_99_percent_Homology", 
                                         "Positive_Epitope_Homology", "Bepipred_Sequence", "detected_epitope", "Bcell_epitope_percentage")


#Change the order of the columns in the dataframe
Post_Filtering_UK1_MHCI_Unique_Peptide_Summary <- Post_Filtering_UK1_MHCI_Unique_Peptide_Summary[, Unique_Peptide_Summary_Column_Order] # leave the row index blank to keep all rows

#Remove the temporary dataframe and variables
remove(temp_df)

####42 Unique MHCI epitopes###

##Create a MHCII table summarizing information per unique peptide

#Summarize the EL Score and Percent Rank for each unique peptide
Post_Filtering_UK1_MHCII_Unique_Peptide_Summary = data.frame()
for (protein in unique(Post_Filtering_UK1_MHCII_Full$Peptide)) {                   #For each unique protein...
  subset_protein = filter(Post_Filtering_UK1_MHCII_Full, Peptide == protein)       #Subset the main dataframe to just that protein
  avg_score_el <- mean(subset_protein$Score_EL)                                    #Get the average Score EL
  avg_percent_rank_el <- mean(subset_protein$X.Rank_EL)                            #Get the average Percent Rank EL
  Post_Filtering_UK1_MHCII_Unique_Peptide_Summary <- rbind(Post_Filtering_UK1_MHCII_Unique_Peptide_Summary, c(protein, avg_score_el, avg_percent_rank_el)) #Save the three elements as a dataframe
}

#Make the column names readable
colnames(Post_Filtering_UK1_MHCII_Unique_Peptide_Summary) <- c("Peptide", "Avg_Score_EL", "Avg_X.Rank_EL")

#Create a dataframe that summarizes all information by Peptide
#This will take unique items in each column and make them into a list
temp_df <- select(Post_Filtering_UK1_MHCII_Full, all_of(Common_MHCI_MHCII_Column_Names)) %>%
  group_by(Peptide) %>%
  # Combine all strings
  summarize_all(funs(toString(unique(.[!is.na(.)])))) %>%
  ungroup()

#Combine the two dataframes while removing repeats and Unique_Id which we will remake
Post_Filtering_UK1_MHCII_Unique_Peptide_Summary <- merge(select(temp_df, -c(Score_EL, X.Rank_EL, Unique_Id)), Post_Filtering_UK1_MHCII_Unique_Peptide_Summary)

#Remake the Unique_Id column with new combined MHC column
Post_Filtering_UK1_MHCII_Unique_Peptide_Summary$Unique_Id <- paste(Post_Filtering_UK1_MHCII_Unique_Peptide_Summary$Identity, "_MHCII_", Post_Filtering_UK1_MHCII_Unique_Peptide_Summary$MHC, "_Pos_", Post_Filtering_UK1_MHCII_Unique_Peptide_Summary$Pos, sep = "")

#Rename the MHC column MHC_Allele
names(Post_Filtering_UK1_MHCII_Unique_Peptide_Summary)[names(Post_Filtering_UK1_MHCII_Unique_Peptide_Summary) == 'MHC'] <- 'MHC_Allele'

#Make a new column called MHC_Type
Post_Filtering_UK1_MHCII_Unique_Peptide_Summary$MHC_Type <- "MHCII"

#Change the order of the columns in the dataframe
Post_Filtering_UK1_MHCII_Unique_Peptide_Summary <- Post_Filtering_UK1_MHCII_Unique_Peptide_Summary[, Unique_Peptide_Summary_Column_Order] # leave the row index blank to keep all rows

#Remove the temporary dataframe and variables
remove(temp_df, avg_percent_rank_el, avg_score_el, subset_protein, protein)

####126 Unique MHCI epitopes###

#Combine with MHCI Table
Post_Filtering_UK1_Unique_Peptide_Summary <- rbind(Post_Filtering_UK1_MHCI_Unique_Peptide_Summary, Post_Filtering_UK1_MHCII_Unique_Peptide_Summary)

#Combine with Protein Annotation Information
Post_Filtering_UK1_Unique_Peptide_Summary <- merge(Post_Filtering_UK1_Unique_Peptide_Summary, Post_Filtering_UK1_Proteins_Annotated_Epitope_Ver)

#Copy this down to Phase II Excel
clipr::write_clip(Post_Filtering_UK1_Unique_Peptide_Summary)

#See Reverse_Vaccinology_UK1_Scripts.txt for next steps for finding Nonoverlapping Epitopes

####Find Nonoverlapping, Unique Epitopes####

#...

#Importing information from Excel...
###41 Unique, Nonoverlapping MHCI with 100% Identity to 99% of All Serovars' Proteomes###
###52 Unique, Nonoverlapping MHCII with 100% Identity to 99% of All Serovars' Proteomes###

#..


###Summarizing Proteins Associated with  Final Epitope Table

##Make Annotated Post Filtering UK1 Proteins version for epitope annotation

#Upload the unique, combined table
Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary <- read.delim("Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary.txt")

#Make no LBL and only LBL versions of the the Nonoverlapping, Unique Epitope Summary
Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary_no_LBL <- filter(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary, LBL_Epitope=="FALSE")
Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary_LBL <- filter(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary, LBL_Epitope=="TRUE")


##Create a verion of the protein annotation dataframe compatible with Epitopes

#Make a copy of Annotated Post Filtering UK1 Proteins
Post_Filtering_UK1_Proteins_Annotated_Epitope_Ver <- Post_Filtering_UK1_Proteins_Annotated

#Change GenBank.Accession to Identity to match Epitope labeling and same columns unique
names(Post_Filtering_UK1_Proteins_Annotated_Epitope_Ver)[names(Post_Filtering_UK1_Proteins_Annotated_Epitope_Ver) == 'GenBank.Accession'] <- 'Identity'
names(Post_Filtering_UK1_Proteins_Annotated_Epitope_Ver)[names(Post_Filtering_UK1_Proteins_Annotated_Epitope_Ver) == 'VaxiJen.Score'] <- 'Protein.VaxiJen.Score'

#Choose columns want to keep
Post_Filtering_UK1_Proteins_Annotated_Epitope_Ver <- select(Post_Filtering_UK1_Proteins_Annotated_Epitope_Ver, c("Identity", "GenBank.Protein.Description", "Similar.Salmonella.Proteins.Descriptions", 
                                                                                                                 "GenBank.Protein.Name", "Similar.Salmonella.Proteins.Protein.Names", 
                                                                                                                 "Vaxign.ML.Score", "Localization", "Localization.Probability", 
                                                                                                                 "Adhesin.Probability", "Trans.membrane.Helices", "Category", 
                                                                                                                 "Protein.VaxiJen.Score", 
                                                                                                                 "Full.Category", "Gene.Name.Essentiality", "Localization.Essentiality", 
                                                                                                                 "Gene.Name.Localization", "RefSeq.Accession", "RefSeq.Protein.Description", 
                                                                                                                 "Ontology.Terms", "GO.Functions"))


##Make a summary of the number of nonoverlapping, unique epitopes associated with each protein

#Aggregate the VaxiJen scores for all epitopes by protein
temp_VaxiJen_Sums_df <- aggregate(as.numeric(VaxiJen_Score) ~ Identity, Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary, sum)
colnames(temp_VaxiJen_Sums_df) <- c("Identity", "Summed_Epitope_VaxiJen_Score")
#53 Proteins with Epitopes

#All Epitope counts by Identity
temp_Total_df <- Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary %>% group_by(Identity) %>% summarise(n_distinct(Peptide))
#53 Proteins with MHCI Epitopes

#Change the column names to something readable
colnames(temp_Total_df) <- c("Identity", "Total_Combined_Unique_Peptides")

#MHCI Peptide counts by Identity
temp_MHCI_df <- filter(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary, MHC_Type == "MHCI") %>% group_by(Identity) %>% summarise(n_distinct(Peptide))
#34 Proteins with MHCI Epitopes

#Change the column names to something readable
colnames(temp_MHCI_df) <- c("Identity", "Combined_Unique_MHCI_Peptides")

#MHCII Peptide counts by Identity
temp_MHCII_df <- filter(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary, MHC_Type == "MHCII") %>% group_by(Identity) %>% summarise(n_distinct(Peptide))
#39 Proteins with MHCI Epitopes

#Change the column names to something readable
colnames(temp_MHCII_df) <- c("Identity", "Combined_Unique_MHCII_Peptides")

#Merge all of the above dataframes into one dataframe
Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary <- merge(temp_VaxiJen_Sums_df, temp_Total_df) 
Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary <- merge(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary, temp_MHCI_df, all.x=TRUE) 
Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary <- merge(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary, temp_MHCII_df, all.x=TRUE) 

#Remove temporary dataframes
remove(temp_MHCI_df, temp_MHCII_df, temp_Total_df, temp_VaxiJen_Sums_df)

#Change any NAs to 0
Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary[is.na(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary)] = 0

#Merge Epitope Counts information with annotation information
Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary <- merge(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary, Post_Filtering_UK1_Proteins_Annotated_Epitope_Ver)
#53 Proteins with either MHCI and/or MHCII epitopes

#Filter the above table by keeping proteins with both MHC Types
Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary_Both_MHC_Types <- rbind(filter(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary, Combined_Unique_MHCI_Peptides > 0 & Combined_Unique_MHCII_Peptides > 0))
#20 Proteins with both MHCI and MHCII epitopes

#Copy down the Protein information to Excel
clipr::write_clip(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Identity_Summary)

####PHASE III CONSTRUCT CREATION####

#Change the working directory
setwd("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Projects_Reverse_Vaccinology/UK1/cleaned_up/Phase_III_Construct_Creation")

###Multiepitope Construct - Most Antigenic Epitopes (9) Localization Agnostic x MHC Type Design

#Find the nine most antigenic epitopes for each MHC Type regardless of location

#MHCI
VaxiJen_MHCI_Localization_Agnostic <- filter(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary_no_LBL, MHC_Type == "MHCI") %>%                                      # Top N highest values by group
  arrange(desc(VaxiJen_Score)) %>% 
  dplyr::slice(1:9)

#MHCI
VaxiJen_MHCII_Localization_Agnostic <- filter(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary_no_LBL, MHC_Type == "MHCII") %>%                                      # Top N highest values by group
  arrange(desc(VaxiJen_Score)) %>% 
  dplyr::slice(1:9)

#Combine the tables
VaxiJen_Localization_Agnostic <- rbind(VaxiJen_MHCI_Localization_Agnostic, VaxiJen_MHCII_Localization_Agnostic, Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary_LBL)

#Chicken Linkers
temp__rand_MHCI_region <- VaxiJen_MHCI_Localization_Agnostic$Peptide %>%
  stringi::stri_rand_shuffle() %>% paste(collapse='AAY')
#Make MHCI region with AAY linker
temp_MHCI_region <- paste(VaxiJen_MHCI_Localization_Agnostic$Peptide, collapse='AAY')

#Make MHCII region with GPGPG linker
temp_MHCII_region <- paste(VaxiJen_MHCII_Localization_Agnostic$Peptide, collapse='GPGPG')

#Make LBL region with KK linker
temp_LBL_region <- paste(Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary_LBL$Peptide, collapse = 'KK')

#Flank epitope region with EAAAK linker
Mugunthan_VaxiJen_Localization_Agnostic_Epitope_Region <- paste('EAAAK', temp_MHCI_region, 'GPGPG', temp_MHCII_region, 'KK', temp_LBL_region, 'EAAAK', sep = "")

#Remove temporary strings
remove(temp_MHCI_region, temp_MHCII_region, temp_LBL_region)

#Determine number of characters
nchar(Mugunthan_VaxiJen_Localization_Agnostic_Epitope_Region)#441


#Copy down both versions of epitope lists and all 4 epitope region sequences to Excel

#Any Localization Constructs
clipr::write_clip(VaxiJen_Localization_Agnostic)
clipr::write_clip(Mugunthan_VaxiJen_Localization_Agnostic_Epitope_Region)

#Make a fasta version of these epitopes
setwd("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Projects_Reverse_Vaccinology/UK1/cleaned_up/Outbreak_Comparisons")
write.fasta(as.list(dput(VaxiJen_Localization_Agnostic$Peptide)), as.list(dput(VaxiJen_Localization_Agnostic$Unique_Id)), "VaxiJen_Localization_Agnostic_Construct_Epitopes.fasta")

####PHASE IV - CONSTRUCT EVALUATION####

####Molecular Dynamics Graphics####

#Load libraries
library(Peptides)
library(tidyverse)
library(scales)
library(ggpubr)

#Save the lengths of construct and TLRs to use later
TLR1 <- 549
TLR2 <- 520
TLR5 <- 441
Loc_Indpdt <- 936

####Localization Independent Docking to TLR1/2 Heteroodimer####

###Pre-Production Run####

##Potential Energy
#Load in the xvg file for Potential Energy
Loc_Indpdt_TLR1_2_potential_E <- readXVG("TLR1_2_Visualizations/potential.xvg")

#Change the Time and Potential to numeric values
Loc_Indpdt_TLR1_2_potential_E$Time <- as.numeric(Loc_Indpdt_TLR1_2_potential_E$Time)
Loc_Indpdt_TLR1_2_potential_E$Potential <- as.numeric(Loc_Indpdt_TLR1_2_potential_E$Potential)

#Plot the data
Loc_Indpdt_TLR1_2_potential_E_line_graph <- ggplot(Loc_Indpdt_TLR1_2_potential_E, aes(x=Time, y=Potential)) +
  geom_line() + 
  #ggtitle("Potential Energy") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(bquote(P[E]~(kJ/mol))) +
  xlab("Time (ps)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR1_2_potential_E_line_graph

tiff('TLR1_2_Visualizations/Loc_Indpdt_TLR1_2_potential_E_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR1_2_potential_E_line_graph
dev.off()

#Record final row of Loc_Indpdt_TLR1_2_potential_E for the manuscript results section

##NVT Temperature
#Load in the xvg file for NVT Temperature
Loc_Indpdt_TLR1_2_nvt_temperature <- readXVG("TLR1_2_Visualizations/nvt_temperature.xvg")

#Change the Time and Temperature to numeric values
Loc_Indpdt_TLR1_2_nvt_temperature$Time <- as.numeric(Loc_Indpdt_TLR1_2_nvt_temperature$Time)
Loc_Indpdt_TLR1_2_nvt_temperature$Temperature <- as.numeric(Loc_Indpdt_TLR1_2_nvt_temperature$Temperature)

#Plot the data
Loc_Indpdt_TLR1_2_nvt_temperature_line_graph <- ggplot(Loc_Indpdt_TLR1_2_nvt_temperature, aes(x=Time, y=Temperature)) +
  geom_line() + 
  #ggtitle("NVT - Temperature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Temperature (K)") +
  xlab("Time (ps)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR1_2_nvt_temperature_line_graph

tiff('TLR1_2_Visualizations/Loc_Indpdt_TLR1_2_nvt_temperature_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR1_2_nvt_temperature_line_graph
dev.off()

#Find the mean throughout the simulation
mean(Loc_Indpdt_TLR1_2_nvt_temperature$Temperature)
#314.619

##NVT RMSD
#Load in the xvg file for NVT RMSD
Loc_Indpdt_TLR1_2_nvt_rmsd <- readXVG("TLR1_2_Visualizations/nvt_rmsd.xvg")

#Change the column names
colnames(Loc_Indpdt_TLR1_2_nvt_rmsd) <- c("Time", "RMSD")

#Change the Time and RMSD to numeric values
Loc_Indpdt_TLR1_2_nvt_rmsd$Time <- as.numeric(Loc_Indpdt_TLR1_2_nvt_rmsd$Time)
Loc_Indpdt_TLR1_2_nvt_rmsd$RMSD <- as.numeric(Loc_Indpdt_TLR1_2_nvt_rmsd$RMSD)

#Plot the data
Loc_Indpdt_TLR1_2_nvt_rmsd_line_graph <- ggplot(Loc_Indpdt_TLR1_2_nvt_rmsd, aes(x=Time, y=RMSD)) +
  geom_line() + 
  #ggtitle("NVT - RMSD") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("RMSD (nm)") +
  xlab("Time (ps)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR1_2_nvt_rmsd_line_graph

tiff('TLR1_2_Visualizations/Loc_Indpdt_TLR1_2_nvt_rmsd_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR1_2_nvt_rmsd_line_graph
dev.off()

#Find the mean throughout the simulation
mean(Loc_Indpdt_TLR1_2_nvt_rmsd$RMSD)
#0.02245825

##NPT Pressure
#Load in the xvg file for NPT Pressure
Loc_Indpdt_TLR1_2_npt_pressure <- readXVG("TLR1_2_Visualizations/npt_pressure.xvg")

#Change the Time and Pressure to numeric values
Loc_Indpdt_TLR1_2_npt_pressure$Time <- as.numeric(Loc_Indpdt_TLR1_2_npt_pressure$Time)
Loc_Indpdt_TLR1_2_npt_pressure$Pressure <- as.numeric(Loc_Indpdt_TLR1_2_npt_pressure$Pressure)

#Plot the data
Loc_Indpdt_TLR1_2_npt_pressure_line_graph <- ggplot(Loc_Indpdt_TLR1_2_npt_pressure, aes(x=Time, y=Pressure)) +
  geom_line() + 
  #ggtitle("NPT - Pressure") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Pressure (bar)") +
  xlab("Time (ps)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR1_2_npt_pressure_line_graph

tiff('TLR1_2_Visualizations/Loc_Indpdt_TLR1_2_npt_pressure_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR1_2_npt_pressure_line_graph
dev.off()

#Find the mean throughout the simulation
mean(Loc_Indpdt_TLR1_2_npt_pressure$Pressure)
#-1.082467

##NPT RMSD
#Load in the xvg file for NPT RMSD
Loc_Indpdt_TLR1_2_npt_rmsd <- readXVG("TLR1_2_Visualizations/npt_rmsd.xvg")

#Change the column names
colnames(Loc_Indpdt_TLR1_2_npt_rmsd) <- c("Time", "RMSD")

#Change the Time and RMSD to numeric values
Loc_Indpdt_TLR1_2_npt_rmsd$Time <- as.numeric(Loc_Indpdt_TLR1_2_npt_rmsd$Time)
Loc_Indpdt_TLR1_2_npt_rmsd$RMSD <- as.numeric(Loc_Indpdt_TLR1_2_npt_rmsd$RMSD)

#Plot the data
Loc_Indpdt_TLR1_2_npt_rmsd_line_graph <- ggplot(Loc_Indpdt_TLR1_2_npt_rmsd, aes(x=Time, y=RMSD)) +
  geom_line() + 
  #ggtitle("NPT - RMSD") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("RMSD (nm)") +
  xlab("Time (ps)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR1_2_npt_rmsd_line_graph

tiff('TLR1_2_Visualizations/Loc_Indpdt_TLR1_2_npt_rmsd_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR1_2_npt_rmsd_line_graph
dev.off()

mean(Loc_Indpdt_TLR1_2_npt_rmsd$RMSD)
#0.0239344


####50 ns Production Run####

##RMSD
#Load in the xvg file for Production Run RMSD
Loc_Indpdt_TLR1_2_50_ns_md_rmsd <- readXVG("TLR1_2_Visualizations/md_50ns_rmsd.xvg")

#Change the column names
colnames(Loc_Indpdt_TLR1_2_50_ns_md_rmsd) <- c("Time", "RMSD")

#Change the Time and RMSD to numeric values
Loc_Indpdt_TLR1_2_50_ns_md_rmsd$Time <- as.numeric(Loc_Indpdt_TLR1_2_50_ns_md_rmsd$Time)
Loc_Indpdt_TLR1_2_50_ns_md_rmsd$RMSD <- as.numeric(Loc_Indpdt_TLR1_2_50_ns_md_rmsd$RMSD)

#Plot the data
Loc_Indpdt_TLR1_2_50_ns_md_rmsd_line_graph <- ggplot(Loc_Indpdt_TLR1_2_50_ns_md_rmsd, aes(x=Time/1000, y=RMSD)) +
  geom_line() + 
  ggtitle("TLR1/2 Heterodimer") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("RMSD (nm)") +
  xlab("Time (ns)") +
  scale_x_continuous(labels = label_comma()) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR1_2_50_ns_md_rmsd_line_graph

tiff('TLR1_2_Visualizations/Loc_Indpdt_TLR1_2_50_ns_md_rmsd_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR1_2_50_ns_md_rmsd_line_graph
dev.off()

##RMSF
#Load in the xvg file for Production Run RMSF
Loc_Indpdt_TLR1_2_50_ns_md_rmsf <- readXVG("TLR1_2_Visualizations/md_50ns_rmsf.xvg")

#Change the column names
colnames(Loc_Indpdt_TLR1_2_50_ns_md_rmsf) <- c("Residue", "RMSF")

#Add a third column based upon which part of the complex each datapoint came by
Loc_Indpdt_TLR1_2_50_ns_md_rmsf$Protein <- c(rep("Construct", Loc_Indpdt), rep("TLR1", 549), rep("TLR2", 520))

#Change the Time and RMSF to numeric values
Loc_Indpdt_TLR1_2_50_ns_md_rmsf$Residue <- as.numeric(Loc_Indpdt_TLR1_2_50_ns_md_rmsf$Residue)
Loc_Indpdt_TLR1_2_50_ns_md_rmsf$RMSF <- as.numeric(Loc_Indpdt_TLR1_2_50_ns_md_rmsf$RMSF)

#Plot the data
Loc_Indpdt_TLR1_2_50_ns_md_rmsf_line_graph <- ggplot(Loc_Indpdt_TLR1_2_50_ns_md_rmsf, aes(x=Residue, y=RMSF,group=Protein)) +
  geom_line(aes(color=Protein)) + 
  #ggtitle("50 ns Production Run - RMSF") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("RMSF (nm)") +
  xlab("Position (bp)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15))
Loc_Indpdt_TLR1_2_50_ns_md_rmsf_line_graph

tiff('TLR1_2_Visualizations/Loc_Indpdt_TLR1_2_50_ns_md_rmsf_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR1_2_50_ns_md_rmsf_line_graph
dev.off()

#Summarize the fluctuations throughout the proteins
Loc_Indpdt_TLR1_2_50_ns_md_rmsf %>%
  group_by(Protein) %>%
  summarise_at(vars(RMSF), list(Min = min, Mean = mean, Sd = sd, Max = max))
# 1 Construct 0.201 0.724 0.360 2.30 
# 2 TLR1      0.120 0.492 0.229 1.20 
# 3 TLR2      0.169 0.458 0.151 0.929

#Get changes in consecutive residues in RMSF
Loc_Indpdt_TLR1_2_50_ns_md_rmsf <- Loc_Indpdt_TLR1_2_50_ns_md_rmsf %>%
  group_by(Protein) %>%
  mutate(Diff = RMSF - lag(RMSF))

#Make a absolute version of the differences
Loc_Indpdt_TLR1_2_50_ns_md_rmsf$Abs_Diff <- abs(Loc_Indpdt_TLR1_2_50_ns_md_rmsf$Diff)

#Summarize changes in consecutive residues in RMSF
Loc_Indpdt_TLR1_2_50_ns_md_rmsf %>%
  na.omit() %>%
  group_by(Protein) %>%
  summarise_at(vars(Abs_Diff), list(Min = min, Mean = mean, Sd = sd, Max = max))
# 1 Construct     0 0.0368 0.0358 0.234 
# 2 TLR1          0 0.0254 0.0184 0.0989
# 3 TLR2          0 0.0235 0.0167 0.172 

##Radius of Gyration
#Load in the xvg file for Production Run Gyrate
Loc_Indpdt_TLR1_2_50_ns_md_gyrate <- readXVG("TLR1_2_Visualizations/md_50ns_gyrate.xvg")

#Change the column names
colnames(Loc_Indpdt_TLR1_2_50_ns_md_gyrate) <- c("Time", "Gyrate")

#Change the Time and Gyrate to numeric values
Loc_Indpdt_TLR1_2_50_ns_md_gyrate$Time <- as.numeric(Loc_Indpdt_TLR1_2_50_ns_md_gyrate$Time)
Loc_Indpdt_TLR1_2_50_ns_md_gyrate$Gyrate <- as.numeric(Loc_Indpdt_TLR1_2_50_ns_md_gyrate$Gyrate)

#Plot the data
Loc_Indpdt_TLR1_2_50_ns_md_gyrate_line_graph <- ggplot(Loc_Indpdt_TLR1_2_50_ns_md_gyrate, aes(x=Time/1000, y=Gyrate)) +
  geom_line() + 
  #ggtitle("50 ns Production Run - Radius of Gyration") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Rg (nm)") +
  xlab("Time (ns)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR1_2_50_ns_md_gyrate_line_graph

tiff('TLR1_2_Visualizations/Loc_Indpdt_TLR1_2_50_ns_md_gyrate_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR1_2_50_ns_md_gyrate_line_graph
dev.off()



####Localization Independent Docking to TLR5 Homodimer####

####Pre-Production Run####

##Potential Energy
#Load in the xvg file for Potential Energy
Loc_Indpdt_TLR5_potential_E <- readXVG("TLR5_Visualizations/potential.xvg")

#Change the Time and Potential to numeric values
Loc_Indpdt_TLR5_potential_E$Time <- as.numeric(Loc_Indpdt_TLR5_potential_E$Time)
Loc_Indpdt_TLR5_potential_E$Potential <- as.numeric(Loc_Indpdt_TLR5_potential_E$Potential)

#Plot the data
Loc_Indpdt_TLR5_potential_E_line_graph <- ggplot(Loc_Indpdt_TLR5_potential_E, aes(x=Time, y=Potential)) +
  geom_line() + 
  #ggtitle("Potential Energy") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(bquote(P[E]~(kJ/mol))) +
  xlab("Time (ps)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR5_potential_E_line_graph

tiff('TLR5_Visualizations/Loc_Indpdt_TLR5_potential_E_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR5_potential_E_line_graph
dev.off()

##NVT Temperature
#Load in the xvg file for NVT Temperature
Loc_Indpdt_TLR5_nvt_temperature <- readXVG("TLR5_Visualizations/nvt_temperature.xvg")

#Change the Time and Temperature to numeric values
Loc_Indpdt_TLR5_nvt_temperature$Time <- as.numeric(Loc_Indpdt_TLR5_nvt_temperature$Time)
Loc_Indpdt_TLR5_nvt_temperature$Temperature <- as.numeric(Loc_Indpdt_TLR5_nvt_temperature$Temperature)

#Plot the data
Loc_Indpdt_TLR5_nvt_temperature_line_graph <- ggplot(Loc_Indpdt_TLR5_nvt_temperature, aes(x=Time, y=Temperature)) +
  geom_line() + 
  #ggtitle("NVT - Temperature") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Temperature (K)") +
  xlab("Time (ps)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR5_nvt_temperature_line_graph

tiff('TLR5_Visualizations/Loc_Indpdt_TLR5_nvt_temperature_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR5_nvt_temperature_line_graph
dev.off()

#Find the mean throughout the simulation
mean(Loc_Indpdt_TLR5_nvt_temperature$Temperature)
#314.6917

##NVT RMSD
#Load in the xvg file for NVT RMSD
Loc_Indpdt_TLR5_nvt_rmsd <- readXVG("TLR5_Visualizations/nvt_rmsd.xvg")

#Change the column names
colnames(Loc_Indpdt_TLR5_nvt_rmsd) <- c("Time", "RMSD")

#Change the Time and RMSD to numeric values
Loc_Indpdt_TLR5_nvt_rmsd$Time <- as.numeric(Loc_Indpdt_TLR5_nvt_rmsd$Time)
Loc_Indpdt_TLR5_nvt_rmsd$RMSD <- as.numeric(Loc_Indpdt_TLR5_nvt_rmsd$RMSD)

#Plot the data
Loc_Indpdt_TLR5_nvt_rmsd_line_graph <- ggplot(Loc_Indpdt_TLR5_nvt_rmsd, aes(x=Time, y=RMSD)) +
  geom_line() + 
  #ggtitle("NVT - RMSD") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("RMSD (nm)") +
  xlab("Time (ps)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR5_nvt_rmsd_line_graph

tiff('TLR5_Visualizations/Loc_Indpdt_TLR5_nvt_rmsd_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR5_nvt_rmsd_line_graph
dev.off()

#Find the mean throughout the simulation
mean(Loc_Indpdt_TLR5_nvt_rmsd$RMSD)
#0.02260239

##NPT Pressure
#Load in the xvg file for NPT Pressure
Loc_Indpdt_TLR5_npt_pressure <- readXVG("TLR5_Visualizations/npt_pressure.xvg")

#Change the Time and Pressure to numeric values
Loc_Indpdt_TLR5_npt_pressure$Time <- as.numeric(Loc_Indpdt_TLR5_npt_pressure$Time)
Loc_Indpdt_TLR5_npt_pressure$Pressure <- as.numeric(Loc_Indpdt_TLR5_npt_pressure$Pressure)

#Plot the data
Loc_Indpdt_TLR5_npt_pressure_line_graph <- ggplot(Loc_Indpdt_TLR5_npt_pressure, aes(x=Time, y=Pressure)) +
  geom_line() + 
  #ggtitle("NPT - Pressure") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Pressure (bar)") +
  xlab("Time (ps)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR5_npt_pressure_line_graph

tiff('TLR5_Visualizations/Loc_Indpdt_TLR5_npt_pressure_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR5_npt_pressure_line_graph
dev.off()

#Find the mean throughout the simulation
mean(Loc_Indpdt_TLR5_npt_pressure$Pressure)
#-3.119037

##NPT RMSD
#Load in the xvg file for NPT RMSD
Loc_Indpdt_TLR5_npt_rmsd <- readXVG("TLR5_Visualizations/npt_rmsd.xvg")

#Change the column names
colnames(Loc_Indpdt_TLR5_npt_rmsd) <- c("Time", "RMSD")

#Change the Time and RMSD to numeric values
Loc_Indpdt_TLR5_npt_rmsd$Time <- as.numeric(Loc_Indpdt_TLR5_npt_rmsd$Time)
Loc_Indpdt_TLR5_npt_rmsd$RMSD <- as.numeric(Loc_Indpdt_TLR5_npt_rmsd$RMSD)

#Plot the data
Loc_Indpdt_TLR5_npt_rmsd_line_graph <- ggplot(Loc_Indpdt_TLR5_npt_rmsd, aes(x=Time, y=RMSD)) +
  geom_line() + 
  #ggtitle("NPT - RMSD") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("RMSD (nm)") +
  xlab("Time (ps)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR5_npt_rmsd_line_graph

tiff('TLR5_Visualizations/Loc_Indpdt_TLR5_npt_rmsd_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR5_npt_rmsd_line_graph
dev.off()

#Find the mean throughout the simulation
mean(Loc_Indpdt_TLR5_npt_rmsd$RMSD)
#0.02388483

####50 ns Production Run####

##RMSD
#Load in the xvg file for Production Run RMSD
Loc_Indpdt_TLR5_50_ns_md_rmsd <- readXVG("TLR5_Visualizations/md_50ns_rmsd_medium.xvg")

#Change the column names
colnames(Loc_Indpdt_TLR5_50_ns_md_rmsd) <- c("Time", "RMSD")

#Change the Time and RMSD to numeric values
Loc_Indpdt_TLR5_50_ns_md_rmsd$Time <- as.numeric(Loc_Indpdt_TLR5_50_ns_md_rmsd$Time)
Loc_Indpdt_TLR5_50_ns_md_rmsd$RMSD <- as.numeric(Loc_Indpdt_TLR5_50_ns_md_rmsd$RMSD)

#Plot the data
Loc_Indpdt_TLR5_50_ns_md_rmsd_line_graph <- ggplot(Loc_Indpdt_TLR5_50_ns_md_rmsd, aes(x=Time/1000, y=RMSD)) +
  geom_line() + 
  ggtitle("TLR5 Homodimer") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("RMSD (nm)") +
  xlab("Time (ns)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR5_50_ns_md_rmsd_line_graph

tiff('TLR5_Visualizations/Loc_Indpdt_TLR5_50_ns_md_rmsd_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR5_50_ns_md_rmsd_line_graph
dev.off()

#Plot the data
Loc_Indpdt_TLR5_50_ns_md_rmsd_line_graph_Chand <- ggplot(Loc_Indpdt_TLR5_50_ns_md_rmsd, aes(x=Time/1000, y=RMSD)) +
  geom_line() + 
  #ggtitle("50 ns Production Run - RMSD") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("RMSD (nm)") +
  xlab("Time (ns)") +
  ylim(0,4.0) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR5_50_ns_md_rmsd_line_graph_Chand

tiff('TLR5_Visualizations/Loc_Indpdt_TLR5_50_ns_md_rmsd_line_graph_Chand.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR5_50_ns_md_rmsd_line_graph_Chand
dev.off()

##RMSF
#Load in the xvg file for Production Run RMSF
Loc_Indpdt_TLR5_50_ns_md_rmsf <- readXVG("TLR5_Visualizations/md_50ns_rmsf_medium.xvg")

#Change the column names
colnames(Loc_Indpdt_TLR5_50_ns_md_rmsf) <- c("Residue", "RMSF")

#Add a third column based upon which part of the complex each datapoint came by
Loc_Indpdt_TLR5_50_ns_md_rmsf$Protein <- c(rep("Construct", Loc_Indpdt), rep("TLR5_A", 441), rep("TLR5_B", 440))

#Change the Time and RMSF to numeric values
Loc_Indpdt_TLR5_50_ns_md_rmsf$Residue <- as.numeric(Loc_Indpdt_TLR5_50_ns_md_rmsf$Residue)
Loc_Indpdt_TLR5_50_ns_md_rmsf$RMSF <- as.numeric(Loc_Indpdt_TLR5_50_ns_md_rmsf$RMSF)

#Plot the data
Loc_Indpdt_TLR5_50_ns_md_rmsf_line_graph <- ggplot(Loc_Indpdt_TLR5_50_ns_md_rmsf, aes(x=Residue, y=RMSF,group=Protein)) +
  geom_line(aes(color=Protein)) + 
  #ggtitle("50 ns Production Run - RMSF") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("RMSF (nm)") +
  xlab("Position (bp)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15))
Loc_Indpdt_TLR5_50_ns_md_rmsf_line_graph

tiff('TLR5_Visualizations/Loc_Indpdt_TLR5_50_ns_md_rmsf_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR5_50_ns_md_rmsf_line_graph
dev.off()

#Summarize the fluctuations throughout the proteins
Loc_Indpdt_TLR5_50_ns_md_rmsf %>%
  group_by(Protein) %>%
  summarise_at(vars(RMSF), list(Min = min, Mean = mean, Sd = sd, Max = max))
# 1 Construct 0.254 0.790 0.247 1.49 
# 2 TLR5_A    0.190 0.375 0.111 0.780
# 3 TLR5_B    0.244 0.565 0.237 1.17 

#Get changes in consecutive residues in RMSF
Loc_Indpdt_TLR5_50_ns_md_rmsf <- Loc_Indpdt_TLR5_50_ns_md_rmsf %>%
  group_by(Protein) %>%
  mutate(Diff = RMSF - lag(RMSF))

#Make a absolute version of the differences
Loc_Indpdt_TLR5_50_ns_md_rmsf$Abs_Diff <- abs(Loc_Indpdt_TLR5_50_ns_md_rmsf$Diff)

#Summarize changes in consecutive residues in RMSF
Loc_Indpdt_TLR5_50_ns_md_rmsf %>%
  na.omit() %>%
  group_by(Protein) %>%
  summarise_at(vars(Abs_Diff), list(Min = min, Mean = mean, Sd = sd, Max = max))
# 1 Construct 0        0.0406 0.0320 0.195 
# 2 TLR5_A    0.000100 0.0272 0.0188 0.108 
# 3 TLR5_B    0.000100 0.0231 0.0165 0.0865

##Radius of Gyration
#Load in the xvg file for Production Run Gyrate
Loc_Indpdt_TLR5_50_ns_md_gyrate <- readXVG("TLR5_Visualizations/md_50ns_gyrate_medium.xvg")

#Change the column names
colnames(Loc_Indpdt_TLR5_50_ns_md_gyrate) <- c("Time", "Gyrate")

#Change the Time and Gyrate to numeric values
Loc_Indpdt_TLR5_50_ns_md_gyrate$Time <- as.numeric(Loc_Indpdt_TLR5_50_ns_md_gyrate$Time)
Loc_Indpdt_TLR5_50_ns_md_gyrate$Gyrate <- as.numeric(Loc_Indpdt_TLR5_50_ns_md_gyrate$Gyrate)

#Plot the data
Loc_Indpdt_TLR5_50_ns_md_gyrate_line_graph <- ggplot(Loc_Indpdt_TLR5_50_ns_md_gyrate, aes(x=Time/1000, y=Gyrate)) +
  geom_line() + 
  #ggtitle("50 ns Production Run - Radius of Gyration") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Rg (nm)") +
  xlab("Time (ns)") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=12), legend.key.size = unit(0.6,"cm"), legend.title=element_text(size=13))
Loc_Indpdt_TLR5_50_ns_md_gyrate_line_graph

tiff('TLR5_Visualizations/Loc_Indpdt_TLR5_50_ns_md_gyrate_line_graph.tiff', units="in", width=16, height=8, res=300)
Loc_Indpdt_TLR5_50_ns_md_gyrate_line_graph
dev.off()

####Make Joint Figures###

#Make a joint figure of pre 50 ns simulation results for TLR1/2 Complex 
Sup_Fig_3_Pre_Prod_TLR1_2_Mol_Dyn <- ggarrange(Loc_Indpdt_TLR1_2_potential_E_line_graph,
                                               Loc_Indpdt_TLR1_2_nvt_rmsd_line_graph,
                                               Loc_Indpdt_TLR1_2_nvt_temperature_line_graph,
                                               Loc_Indpdt_TLR1_2_npt_rmsd_line_graph,
                                               Loc_Indpdt_TLR1_2_npt_pressure_line_graph,
                                               labels = c("A", "B", "C", "D", "E"),
                                               ncol = 2, nrow = 3)
Sup_Fig_3_Pre_Prod_TLR1_2_Mol_Dyn

tiff('Supplementary Figure 3 - Pre-Production TLR1-2 Molecular Dynamics - GROMACS.tiff', units="in", width=8, height=8, res=300)
Sup_Fig_3_Pre_Prod_TLR1_2_Mol_Dyn
dev.off()

#Make a joint figure of pre 50 ns simulation results for TLR5 Complex 
Sup_Fig_4_Pre_Prod_TLR5_Mol_Dyn <- ggarrange(Loc_Indpdt_TLR5_potential_E_line_graph,
                                             Loc_Indpdt_TLR5_nvt_rmsd_line_graph,
                                             Loc_Indpdt_TLR5_nvt_temperature_line_graph,
                                             Loc_Indpdt_TLR5_npt_rmsd_line_graph,
                                             Loc_Indpdt_TLR5_npt_pressure_line_graph,
                                             labels = c("A", "B", "C", "D", "E"),
                                             ncol = 2, nrow = 3)
Sup_Fig_4_Pre_Prod_TLR5_Mol_Dyn

tiff('Supplementary Figure 4 - Pre-Production TLR5 Molecular Dynamics - GROMACS.tiff', units="in", width=8, height=8, res=300)
Sup_Fig_4_Pre_Prod_TLR5_Mol_Dyn
dev.off()


#Make a joint figure of 50 ns simulations results 
Figure_5_Mol_Dyn <- ggarrange(Loc_Indpdt_TLR1_2_50_ns_md_rmsd_line_graph,
                              Loc_Indpdt_TLR5_50_ns_md_rmsd_line_graph,
                              Loc_Indpdt_TLR1_2_50_ns_md_rmsf_line_graph,
                              Loc_Indpdt_TLR5_50_ns_md_rmsf_line_graph,
                              Loc_Indpdt_TLR1_2_50_ns_md_gyrate_line_graph,
                              Loc_Indpdt_TLR5_50_ns_md_gyrate_line_graph,
                              labels = c("A", "B", "C", "D", "E", "F"),
                              ncol = 2, nrow = 3)
Figure_5_Mol_Dyn



tiff('Figure 5 - Molecular Dynamics - GROMACS.tiff', units="in", width=8, height=8, res=300)
Figure_5_Mol_Dyn
dev.off()

####Epitope Homology to PulseNet-Outbreak Isolates

#Load libraries
library(tidyverse)
library(lubridate)
library(foreach)
library(parallel)
library(doParallel)

####PulseNet 030124 Metadata####

#Set the working directory
setwd("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Bioinformatics_Analysis/PulseNet_030124")

#Prep for running in parallel

#Detect number of available cores and create cluster
cl <- parallel::makeCluster(6)

# Activate cluster for foreach library
doParallel::registerDoParallel(cl)

#Upload the PulseNet isolate information
#One isolate failed year parsing, change it from 2021-10-10- to 10/10/2021
PulseNet_Isolate_Information_V2 <- read.csv("PulseNet_Isolate_Information_V2.csv", na.strings=c(""," ","NA"))

#Upload the PulseNet outbreak information
#Some Outbreaks has multiple values, edit as appropriate
PulseNet_Outbreak_Information <- read.csv("PulseNet_Outbreak_Info.csv", na.strings=c(""," ","NA"))
colnames(PulseNet_Outbreak_Information) <- c("Outbreak", "Outbreak_Source")

#Merge the isolate and outbreak information
PulseNet_Full_V2 <- merge(PulseNet_Isolate_Information_V2, PulseNet_Outbreak_Information, all.x = TRUE)

#Create a collection year column
PulseNet_Full_V2$Collection.year <- year(parse_date_time(PulseNet_Full_V2$IsolatDate, orders = c("m/d/y")))

#Change the order of columns for easy reading
PulseNet_Full_V2<- subset(PulseNet_Full_V2, select =  c("WGS_id", "IsolatDate", "Collection.year", "NCBI_ACCESSION", "PatientAgeDays", 
                                                        "PatientAgeMonths", "PatientAgeYears", "PatientSex", "PulseNet_UploadDate", 
                                                        "ReceivedDate", "REP_code", "Serotype_wgs", "SourceCity", "SourceCountry", 
                                                        "SourceCounty", "SourceSite", "SourceState", "SourceType", "Type.Details",  "SRR_id", "Outbreak",
                                                        "Outbreak_Source"))

#Create a table of the distribution of isolation sources
PulseNet_Full_V2_by_isolate_source <- PulseNet_Full_V2 %>%
  dplyr::count(SourceSite, SourceType, Type.Details)

##Add Summarized Isolation Source Column

#Combine the Isolation.source and Host information to get a column to search for patterns for Isolation Source Summary column
PulseNet_Full_V2 <- unite(PulseNet_Full_V2, "SourceSite_SourceType_TypeDetails", c(SourceSite, SourceType, Type.Details), remove=FALSE, sep = "_")

#Make lists of terms associated with organisms of interest
Turkey_ids <- c("Turkey", "turkey", "Meleagris gallopavo", "meleagris gallopavo")

Chicken_ids <- c("Chicken", "chicken", "Gallus gallus", "Gallus Gallus", "CHICKEN", "chick", "Chick", "CHICK PADS", "broiler", "Broiler", "Layers", "layers")

Cow_ids <- c("Cow", "cow", "Cattle", "cattle", "Bos taurus", "beef", "Beef", "bovine", "Bovine", "Calf", "calf")

Swine_ids <- c("Swine", "swine", "Sus scrofa", "porcine", "Porcine", "pork", "Pork", "PORK", "Hogs", "pig", "Pig", "Sow", "sow")

Human_ids <- c("Human", "human", "HUMAN", "Homo sapiens", "homo sapiens")

Egg_ids <- c("egg", "EGG", "Egg")

#Note: If an entry did not contain any of the above, but said poultry it was put into the Other Source/Host since poultry could mean turkey or chicken

#Use a for loop to look at each row and summarize Isolation/host information
PulseNet_Full_V2_Summ_Isolation_Source <- foreach(row=1:nrow(PulseNet_Full_V2), .combine = rbind, .packages=c('tidyverse')) %dopar% {                                         #For each row in the metadata file...
  WGS_id <- PulseNet_Full_V2[row, "WGS_id"]                         #Save the WGS_id
  SourceSite_SourceType_TypeDetails <- PulseNet_Full_V2[row, "SourceSite_SourceType_TypeDetails"]
  if (grepl(paste(Turkey_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE){                                      #If the Run Number is in the pESI positive list
    Summarized.isolation.source = "Turkey-Associated"                                                                 #Label the pESI_Status as Present
  }
  else if (grepl(paste(Chicken_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Chicken-Associated"
  }
  else if (grepl(paste(Cow_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Cow-Associated"
  }
  else if (grepl(paste(Swine_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Swine-Associated"
  }
  else if (grepl(paste(Human_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Human-Associated"
  }
  else if (grepl("NA_NA_NA", SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Not Indicated"
  }
  else{                                                                                    #If the Run Number is not in the pESI positive list
    Summarized.isolation.source = "Other Source/Host"                                                             #Label the pESI_Status as Not Present
  }
  data.frame(WGS_id, SourceSite_SourceType_TypeDetails, Summarized.isolation.source)  #Save the three elements per row as a dataframe
}

PulseNet_Full_V2_Summ_Isolation_Source <- foreach(row=1:nrow(PulseNet_Full_V2), .combine = rbind, .packages=c('tidyverse')) %dopar% {                                         #For each row in the metadata file...
  WGS_id <- PulseNet_Full_V2[row, "WGS_id"]                         #Save the WGS_id
  SourceSite_SourceType_TypeDetails <- PulseNet_Full_V2[row, "SourceSite_SourceType_TypeDetails"]
  if (grepl(paste(Egg_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE){                                      #If the Run Number is in the pESI positive list
    Summarized.isolation.source = "Egg-Associated"                                                                 #Label the pESI_Status as Present
  }
  else if (grepl(paste(Turkey_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Turkey-Associated"
  }
  else if (grepl(paste(Chicken_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Chicken-Associated"
  }
  else if (grepl(paste(Cow_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Cow/Beef-Associated"
  }
  else if (grepl(paste(Swine_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Swine/Pork-Associated"
  }
  else if (grepl(paste(Human_ids, collapse = "|"), SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Human-Associated"
  }
  else if (grepl("clinical", SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Other/Unspecified Clinical"
  }
  else if (grepl("NA_NA_NA", SourceSite_SourceType_TypeDetails) == TRUE) {
    Summarized.isolation.source = "Not Indicated"
  }
  else{                                                                                    #If the Run Number is not in the pESI positive list
    Summarized.isolation.source = "Other Source/Host"                                                             #Label the pESI_Status as Not Present
  }
  data.frame(WGS_id, SourceSite_SourceType_TypeDetails, Summarized.isolation.source)  #Save the three elements per row as a dataframe
}


#After draft is complete check for ids to potentially summarize
temp_PulseNet_Full_V2_Summ_Isolation_Source_Other <- unique(filter(PulseNet_Full_V2_Summ_Isolation_Source, Summarized.isolation.source=="Other Source/Host")$SourceSite_SourceType_TypeDetails)
temp_PulseNet_Full_V2_Summ_Isolation_Source_Other
#314 unique labels in Other Source/Host

#Draft search function
dput(grep("Pig", temp_PulseNet_Full_V2_Summ_Isolation_Source_Other, value = TRUE))

#Merge with the full metadata dataframe
PulseNet_Full_V2_summ <- merge(PulseNet_Full_V2, PulseNet_Full_V2_Summ_Isolation_Source)

#Check how well the main metadata columns match one another

table(PulseNet_Full_V2_summ$Summarized.isolation.source == PulseNet_Full_V2_summ$Summarized.isolation.source)
# FALSE  TRUE 
# 3509 30255

table(filter(PulseNet_Full_summ, Serotype_wgs == "Infantis")$Summarized.isolation.source == filter(PulseNet_Full_V2_summ, Serotype_wgs == "Infantis")$Summarized.isolation.source)
# FALSE  TRUE 
# 81  1875

table(filter(PulseNet_Full_summ, Serotype_wgs == "Senftenberg")$Summarized.isolation.source == filter(PulseNet_Full_V2_summ, Serotype_wgs == "Senftenberg")$Summarized.isolation.source)
# FALSE  TRUE 
# 4   127

#Add a column combining location information
PulseNet_Full_V2_summ <- unite(PulseNet_Full_V2_summ, "Location", c(SourceCountry, SourceState, SourceCounty,  SourceCity), remove=FALSE, sep = ":")

#Get a list of all SourceState and double check that they are all actual USA States or Territories and then save it
dput(sort(unique(PulseNet_Full_V2_summ$SourceState)))
USA_States_Territories <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DC", "DE", "FL", 
                            "Ga", "GA", "Guam", "HI", "IA", "ID", "IL", "IN", "KS", "KS, TX, CO", 
                            "KY", "LA", "MA", "MD", "ME", "MI", "Mn", "MN", "MO", "MS", "MT", 
                            "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "NYC", "OH", 
                            "OK", "OR", "PA", "PR", "RI", "SC", "SD", "TN", "TX", "UT", "VA", 
                            "VT", "WA", "WI", "WV", "WY")


##Add Summarized Isolation Location Column
#Based upon the following wikipedia page
#https://simple.wikipedia.org/wiki/List_of_countries_by_continents

#Upload country names dataframe
Country_Names <- read.delim("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Bioinformatics_Analysis/Country_Names.txt")

#Save each column as a list
Africa_ids_wiki <- na.omit(Country_Names$Africa)

Asia_ids_wiki <- na.omit(Country_Names$Asia)

Europe_ids_wiki <- na.omit(Country_Names$Europe)

Oceania_ids_wiki <- na.omit(Country_Names$Oceania)

South_America_ids_wiki <- na.omit(Country_Names$South.America)

Central_America_ids_wiki <- na.omit(Country_Names$Other.North.America)

#Make lists of terms associated with organisms of interest
United_States_ids <- c("USA", "United States", "Puerto Rico")

United_Kingdom_ids <- c("United Kingdom", "Ireland", "England")

Europe_ids <- c(Europe_ids_wiki, "Czech Republic")

Asia_ids <- c(Asia_ids_wiki, "Korea", "West Bank", "Viet Nam", "Hong Kong")

Oceania_ids <- c(Oceania_ids_wiki)

Central_America_Caribbean <- c(Central_America_ids_wiki, "Bahamas", "Cayman Islands")

Canada_ids <- c("Canada")

South_America_ids <- c(South_America_ids_wiki)

Africa_ids <- c(Africa_ids_wiki)

#Use a for loop to look at each row and summarize Isolation/host information
PulseNet_Full_V2_summ_Summ_Location <- foreach(row=1:nrow(PulseNet_Full_V2_summ), .combine = rbind, .packages=c('tidyverse')) %dopar% {                                     #For each row in the metadata file...
  WGS_id <- PulseNet_Full_V2_summ[row, "WGS_id"]                         #Save the WGS_id
  Location <- PulseNet_Full_V2_summ[row, "Location"]
  SourceState <- PulseNet_Full_V2_summ[row, "SourceState"]
  if (grepl(paste(United_States_ids, collapse = "|"), Location) == TRUE){                                      #If the Run Number is in the pESI positive list
    Summarized.location = "North America: USA"                                                                 #Label the pESI_Status as Present
  }
  else if (grepl(paste(Europe_ids, collapse = "|"), Location) == TRUE) {
    Summarized.location = "Europe: Other"
  }
  else if (grepl(paste(United_Kingdom_ids, collapse = "|"), Location) == TRUE) {
    Summarized.location = "Europe: United Kingdom"
  }
  else if (grepl(paste(Oceania_ids, collapse = "|"), Location) == TRUE) {
    Summarized.location = "Oceania"
  }
  else if (grepl(paste(Asia_ids, collapse = "|"), Location) == TRUE) {
    Summarized.location = "Asia"
  }
  else if (grepl(paste(Central_America_Caribbean, collapse = "|"), Location) == TRUE) {
    Summarized.location = "North America: Other"
  }
  else if (grepl(paste(Canada_ids, collapse = "|"), Location) == TRUE) {
    Summarized.location = "North America: Canada"
  }
  else if (grepl(paste(South_America_ids, collapse = "|"), Location) == TRUE) {
    Summarized.location = "South America"
  }
  else if (grepl(paste(Africa_ids, collapse = "|"), Location) == TRUE) {
    Summarized.location = "Africa"
  }
  else if (grepl(paste(USA_States_Territories, collapse = "|"), SourceState) == TRUE) {
    Summarized.location = "North America: USA"
  }
  else if (grepl(paste(Europe_ids, collapse = "|"), Location) == TRUE) {
    Summarized.location = "Europe: Other"
  }
  else{                                                                                    #If the Run Number is not in the pESI positive list
    Summarized.location = "Not Indicated"                                                             #Label the pESI_Status as Not Present
  }
  data.frame(WGS_id, Location, Summarized.location)  #Save the four elements per row as a dataframe
}


#After draft is complete check for ids to potentially summarize
temp_PulseNet_Full_V2_summ_Summ_Location_Other <- unique(filter(PulseNet_Full_V2_summ_Summ_Location, Summarized.location=="Not Indicated")$Location)
temp_PulseNet_Full_V2_summ_Summ_Location_Other

#Merge with the full metadata dataframe
PulseNet_Full_V2_summ <- merge(PulseNet_Full_V2_summ, PulseNet_Full_V2_summ_Summ_Location)

##Create a North America: USA Location Census Geographic Regions dataframe to combine with the larger metadata file

#Make a USA focused dataframe
PulseNet_Full_V2_summ_USA <- dplyr::filter(PulseNet_Full_V2_summ, Summarized.location=="North America: USA")
#33323

#Make lists of terms associated with organisms of interest
CGR_Northeast_States <- c("CT", "ME", "MA", "NH", "RI", "VT",
                          "NJ", "NY", "PA")

CGR_Midwest_States <- c("IL", "IN", "MI", "OH", "WI",
                        "IA", "KS", "MN", "Mn", "MO", "NE", "ND", "SD")

CGR_South_States <- c("DE", "DC", "FL", "GA", "Ga", "MD", "NC", "SC", "VA", "WV", "PR",
                      "AL", "KY", "MS", "TN",
                      "AR", "LA", "OK", "TX")

CGR_West_States <- c("AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY",
                     "AK", "CA", "HI", "OR", "WA", "Guam")

#Use a for loop to look at each row and summarize Isolation/host information
PulseNet_Full_V2_summ_USA_Region <- foreach(row=1:nrow(PulseNet_Full_V2_summ_USA), .combine = rbind, .packages=c('tidyverse')) %dopar% {                                         #For each row in the metadata file...
  WGS_id <- PulseNet_Full_V2_summ_USA[row, "WGS_id"]                         #Save the WGS_id
  Location <- PulseNet_Full_V2_summ_USA[row, "Location"]
  SourceState <- PulseNet_Full_V2_summ_USA[row, "SourceState"]
  if (grepl(paste(CGR_Northeast_States, collapse = "|"), SourceState) == TRUE){                                      #If the Run Number is in the pESI positive list
    USA.region.summarized.location = "USA: Northeast Region"                                                                 #Label the pESI_Status as Present
  }
  else if (grepl(paste(CGR_Midwest_States, collapse = "|"), SourceState) == TRUE) {
    USA.region.summarized.location = "USA: Midwest Region"
  }
  else if (grepl(paste(CGR_South_States, collapse = "|"), SourceState) == TRUE) {
    USA.region.summarized.location = "USA: South Region"
  }
  else if (grepl(paste(CGR_West_States, collapse = "|"), SourceState) == TRUE) {
    USA.region.summarized.location = "USA: West Region"
  }
  else{                                                                                    #If the Run Number is not in the pESI positive list
    USA.region.summarized.location = "USA: Undefined Region"                                                             #Label the pESI_Status as Not Present
  }
  data.frame(WGS_id, Location, SourceState, USA.region.summarized.location)  #Save the four elements per row as a dataframe
}

#After draft is complete check for ids to potentially summarize
unique(filter(PulseNet_Full_V2_summ_USA_Region, USA.region.summarized.location=="USA: Undefined Region")$Location)
#17

#Merge with the full metadata dataframe
PulseNet_Full_V2_summ_USA <- merge(PulseNet_Full_V2_summ_USA, PulseNet_Full_V2_summ_USA_Region, by = c("WGS_id", "Location", "SourceState"))

##Create a North America: USA Location Census Geographic Regions dataframe to combine with the larger metadata file
#Make lists of terms associated with organisms of interest
CGD_New_England_States <- c("CT", "ME", "MA", "NH", "RI", "VT") 
CGD_Middle_Atlantic_States <- c("NJ", "NY", "PA")

CGD_East_North_Central_States <- c("IL", "IN", "MI", "OH", "WI")
CGD_West_North_Central_States <- c("IA", "KS", "MN", "Mn", "MO", "NE", "ND", "SD")

CGD_South_Atlantic_States <- c("DE", "DC", "FL", "GA", "Ga", "MD", "NC", "SC", "VA", "WV", "PR")
CGD_East_South_Central_States <- c("AL", "KY", "MS", "TN")
CGD_West_South_Central_States <- c("AR", "LA", "OK", "TX")

CGD_Mountain_States <- c("AZ", "CO", "ID", "MT", "NV", "NM", "UT", "WY")
CGD_Pacific_States <- c("AK", "CA", "HI", "OR", "WA", "Guam")

#Use a for loop to look at each row and summarize Isolation/host information
PulseNet_Full_V2_summ_USA_Division <- foreach(row=1:nrow(PulseNet_Full_V2_summ_USA), .combine = rbind, .packages=c('tidyverse')) %dopar% {                                         #For each row in the metadata file...
  WGS_id <- PulseNet_Full_V2_summ_USA[row, "WGS_id"]                         #Save the WGS_id
  Location <- PulseNet_Full_V2_summ_USA[row, "Location"]
  SourceState <- PulseNet_Full_V2_summ_USA[row, "SourceState"]
  if (grepl(paste(CGD_New_England_States, collapse = "|"), SourceState) == TRUE){                                      #If the Run Number is in the pESI positive list
    USA.division.summarized.location = "USA: New England Division"                                                                 #Label the pESI_Status as Present
  }
  else if (grepl(paste(CGD_Middle_Atlantic_States, collapse = "|"), SourceState) == TRUE) {
    USA.division.summarized.location = "USA: Middle Atlantic Division"
  }
  else if (grepl(paste(CGD_East_North_Central_States, collapse = "|"), SourceState) == TRUE) {
    USA.division.summarized.location = "USA: East North Central Division"
  }
  else if (grepl(paste(CGD_West_North_Central_States, collapse = "|"), SourceState) == TRUE) {
    USA.division.summarized.location = "USA: West North Central Division"
  }
  else if (grepl(paste(CGD_South_Atlantic_States, collapse = "|"), SourceState) == TRUE) {
    USA.division.summarized.location = "USA: South Atlantic Division"
  }
  else if (grepl(paste(CGD_East_South_Central_States, collapse = "|"), SourceState) == TRUE) {
    USA.division.summarized.location = "USA: East South Central Division"
  }
  else if (grepl(paste(CGD_West_South_Central_States, collapse = "|"), SourceState) == TRUE) {
    USA.division.summarized.location = "USA: West South Central Division"
  }
  else if (grepl(paste(CGD_Mountain_States, collapse = "|"), SourceState) == TRUE) {
    USA.division.summarized.location = "USA: Mountain Division"
  }
  else if (grepl(paste(CGD_Pacific_States, collapse = "|"), SourceState) == TRUE) {
    USA.division.summarized.location = "USA: Pacific Division"
  }
  else{                                                                                    #If the Run Number is not in the pESI positive list
    USA.division.summarized.location = "USA: Undefined Division"                                                             #Label the pESI_Status as Not Present
  }
  data.frame(WGS_id, Location, SourceState, USA.division.summarized.location)  #Save the four elements per row as a dataframe
}


#After draft is complete check for ids to potentially summarize
unique(filter(PulseNet_Full_V2_summ_USA_Division, USA.division.summarized.location=="USA: Undefined Division")$Location)

#Merge with the full metadata dataframe
PulseNet_Full_V2_summ_USA <- merge(PulseNet_Full_V2_summ_USA, PulseNet_Full_V2_summ_USA_Division)

#Make the Collection.year column numeric
PulseNet_Full_V2_summ$Collection.year <- as.numeric(PulseNet_Full_V2_summ$Collection.year)
PulseNet_Full_V2_summ_USA$Collection.year <- as.numeric(PulseNet_Full_V2_summ_USA$Collection.year)

#Reorder USA Region and Division columns in USA-focused dataframe
PulseNet_Full_V2_summ_USA$USA.region.summarized.location <- factor(PulseNet_Full_V2_summ_USA$USA.region.summarized.location, levels=c("USA: Northeast Region", "USA: Midwest Region", "USA: South Region", "USA: West Region", "USA: Undefined Region"))
PulseNet_Full_V2_summ_USA$USA.division.summarized.location <- factor(PulseNet_Full_V2_summ_USA$USA.division.summarized.location, levels=c("USA: New England Division", "USA: Middle Atlantic Division", "USA: East North Central Division", "USA: West North Central Division", "USA: South Atlantic Division", "USA: East South Central Division", "USA: West South Central Division", "USA: Mountain Division", "USA: Pacific Division", "USA: Undefined Division"))

#Merge wtih full metadata dataframe
PulseNet_Full_V2_summ <- merge(PulseNet_Full_V2_summ, PulseNet_Full_V2_summ_USA, all = TRUE)


#Create a table of the distribution of whole genome defined serovars
PulseNet_Full_V2_summ_by_targets <- PulseNet_Full_V2_summ %>%
  dplyr::count(Serotype_wgs, Summarized.isolation.source, Collection.year, Summarized.location) %>%
  group_by(Serotype_wgs, Summarized.isolation.source, Collection.year, Summarized.location)

#Create a table of the distribution of whole genome defined serovars by source and total
PulseNet_Full_V2_by_serovar <- PulseNet_Full_V2_summ %>%
  dplyr::count(Serotype_wgs, Summarized.isolation.source) %>%
  group_by(Serotype_wgs)%>%
  tidyr::pivot_wider(names_from = Summarized.isolation.source, values_from = n, 
                     values_fill = list(n = 0))

PulseNet_Full_V2_by_serovar$Total <- rowSums(PulseNet_Full_V2_by_serovar[,2:length(PulseNet_Full_V2_by_serovar)])


Serogroups <- read.csv("Salmonella-serotype_serogroup_antigen_table-WHO_2007.csv")

PulseNet_Full_V2_by_serovar <- merge(PulseNet_Full_V2_by_serovar, Serogroups, by.x = "Serotype_wgs", by.y = "Serovar", all.x = TRUE)

#Create a table of the distribution of whole genome defined serovars by REP
PulseNet_Full_V2_by_REP_code <- filter(PulseNet_Full_V2_summ, !is.na(REP_code)) %>%
  dplyr::count(Serotype_wgs, REP_code) 

####Salmonella Multiepitope Vaccine Outbreak Comparisons####

#Change the directory
setwd("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Projects_Reverse_Vaccinology/UK1/cleaned_up/Outbreak_Comparisons")

##Enteritidis

#Make a list of all Senftenberg associated SRRs
PulseNet_Enteritidis_SRRs <- unique(filter(PulseNet_Full_V2_summ, Serotype_wgs =="Enteritidis")$SRR_id)
#6245

#Save list as a txt file and send to HPC
write_delim(as.data.frame(PulseNet_Enteritidis_SRRs), "PulseNet_Enteritidis_SRRs.txt", col_names = FALSE)

##Newport & Typhimurium 

#Make a list of all Newport and Typhimurium associated SRRs
PulseNet_Newport_Typhimurium_SRRs <- unique(filter(PulseNet_Full_V2_summ, Serotype_wgs %in% c("Newport", "Typhimurium"))$SRR_id)
#7767

#Save list as a txt file and send to HPC
write_delim(as.data.frame(PulseNet_Newport_Typhimurium_SRRs), "PulseNet_Newport_Typhimurium_SRRs.txt", col_names = FALSE)

## Less than 2000 but Greater than 1000 = Hadar, Infantis (not in above), Braenderup, Oranienburg, I 4:i:-

#Make a list of all Serotypes besides Enteritidis, Newport, and Typhimurium with >1000 isolates associated SRRs
PulseNet_other_GT_1000_SRRs <- unique(filter(PulseNet_Full_V2_summ, Serotype_wgs %in% c("Hadar", "Infantis", "Braenderup", "Oranienburg", "I 4:i:-"))$SRR_id)
#7935

#Make list of Infantis isolates that have already been downloaded due to Senftenberg study (n=1008)
PulseNet_Infantis_19_23_Human_Turkey_SRRs <- intersect(PulseNet_other_GT_1000_SRRs, PulseNet_V2_Infantis_USA_2019_2023_021524_Human_Turkey$SRR_id)
#1008

#Save list as a txt file and send to HPC
write_delim(as.data.frame(PulseNet_Infantis_19_23_Human_Turkey_SRRs), "PulseNet_Infantis_19_23_Human_Turkey_SRRs.txt", col_names = FALSE)

#Remove Infantis isolates that have already been downloaded due to Senftenberg study (n=1008)
PulseNet_other_GT_1000_SRRs <- setdiff(PulseNet_other_GT_1000_SRRs, PulseNet_V2_Infantis_USA_2019_2023_021524_Human_Turkey$SRR_id)
#6927

#Save list as a txt file and send to HPC
write_delim(as.data.frame(PulseNet_other_GT_1000_SRRs), "PulseNet_other_GT_1000_SRRs.txt", col_names = FALSE)



## Less than 1000 but Greater than 250 = Berta, Heidelberg, Mbandaka, Montevideo, Muenchen, Paratyphi B var. L(+) tartrate+, Poona, Reading, Saintpaul, Sundsvall, Thompson, Uganda

#Make a list of all Newport and Typhimurium associated SRRs
PulseNet_LT_1000_GT_250_SRRs <- unique(filter(PulseNet_Full_V2_summ, Serotype_wgs %in% dput(filter(PulseNet_Full_V2_by_serovar, Total > 250 & Total < 1000)$Serotype_wgs))$SRR_id)
#5763

#Save list as a txt file and send to HPC
write_delim(as.data.frame(PulseNet_LT_1000_GT_250_SRRs), "PulseNet_LT_1000_GT_250_SRRs.txt", col_names = FALSE)

## Less than 250 (Serovars=122)

#Make a list of all Newport and Typhimurium associated SRRs
PulseNet_LT_250_SRRs <- unique(filter(PulseNet_Full_V2_summ, Serotype_wgs %in% dput(filter(PulseNet_Full_V2_by_serovar, Total < 250)$Serotype_wgs))$SRR_id)
#5039

#Save list as a txt file and send to HPC
write_delim(as.data.frame(PulseNet_LT_250_SRRs), "PulseNet_LT_250_SRRs.txt", col_names = FALSE)

#Make a list of all isolates
PulseNet_All_SRRs <- unique(PulseNet_Full_V2_summ$SRR_id)
#33763

write_delim(as.data.frame(PulseNet_All_SRRs), "PulseNet_All_SRRs.txt", col_names = FALSE)

#Change working directory back
setwd("C:/Users/David.Bradshaw/OneDrive - USDA/Documents/Bioinformatics_Analysis/PulseNet_030124")

##CheckM

#Upload CheckM results
PulseNet_All_CheckM <- read.delim("PulseNet_checkm_results_edited.txt", header = TRUE)

#Make a new column that indicates if an assembly passes the filters
#Completeness >= 95% & Contamination <= 5%
PulseNet_All_CheckM$CheckM_Pass <- ifelse(PulseNet_All_CheckM$Completeness>=95 & PulseNet_All_CheckM$Contamination<=5, TRUE, FALSE)

#Get the distribution
table(PulseNet_All_CheckM$CheckM_Pass)
# FALSE  TRUE 
# 82    33672 
#Lose 82/33754 isolates (0.24%)

#Make both versions of CheckM passing or not passing
PulseNet_All_CheckM_Pass <- PulseNet_All_CheckM_No_Pass <- filter(PulseNet_All_CheckM, CheckM_Pass==TRUE)
PulseNet_All_CheckM_No_Pass <- filter(PulseNet_All_CheckM, CheckM_Pass==FALSE)


#Upload the proteomes epitope homology Step 1 table
PulseNet_proteomes_epitope_homology_Step_1 <- read.delim("PulseNet_proteomes_epitope_homology.txt", header = TRUE, sep = " ")

#Check number of isolates
length(unique(PulseNet_proteomes_epitope_homology_Step_1$Proteome_Name))
#33754

#Remove isolates that do not pass the CheckM filter
PulseNet_proteomes_epitope_homology_Step_1_QCed <- filter(PulseNet_proteomes_epitope_homology_Step_1, Proteome_Name %in% filter(PulseNet_All_CheckM, CheckM_Pass==TRUE)$WGS_id)

#Check that it worked
length(unique(PulseNet_proteomes_epitope_homology_Step_1_QCed$Proteome_Name))
#33672

#Send it to HPC
write_delim(PulseNet_proteomes_epitope_homology_Step_1_QCed, "PulseNet_proteomes_epitope_homology_Step_1_QCed.txt")

##Post epitope - proteome blasting summary

#Upload results
PulseNet_proteomes_epitope_homology_summary <- read.delim("PulseNet_proteomes_epitope_homology_summary.txt", header = TRUE)

#Get distibrution of epitopes with 100% identity and coverage (lengths of 9 or 15 for MHCI or MCHII respectively)
#with homology to 99% of PulseNet Proteomes
table(PulseNet_proteomes_epitope_homology_summary$Full_0_GT_99)
# FALSE  TRUE 
# 9    19

#Get distibrution of epitopes with 100% identity and coverage (lengths of 9 or 15 for MHCI or MCHII respectively)
#with homology to 90% of PulseNet Proteomes
table(PulseNet_proteomes_epitope_homology_summary$Full_0_GT_90)
# FALSE  TRUE 
# 3    25

#Get distibrution of epitopes with >= 88% identity and coverage (lengths of 9/8 or 15/14 for MHCI or MCHII respectively)
#with homology to 90% of PulseNet Proteomes
table(PulseNet_proteomes_epitope_homology_summary$Partial_1_GT_90)
# FALSE  TRUE 
# 1    27

####Overlap with Positive Homology Isolates####

#Upload list of all Enteritidis, Typhimurium, Infantis, Kentucky, Hadar, and Uganda isolates from NCBI PDD (by Computed Types)
PDD_Target_Serovars_011025_Metadata <- read.delim("Target_Serovars_011025.tsv", na.strings=c(""," ","NA"))
#266000

#Separate out Computed.types to create a IB_SeqSero2.serotype column
PDD_Target_Serovars_011025_Metadata <- separate(PDD_Target_Serovars_011025_Metadata, Computed.types, sep = "=", into = c("Extra1", "Extra2", "IB_SeqSero2.serotype"), remove = FALSE)

#Get distribution
table(PDD_Target_Serovars_011025_Metadata$IB_SeqSero2.serotype, useNA = "always")
# Enteritidis       Hadar    Infantis    Kentucky Typhimurium      Uganda        <NA> 
#   122209        5777       34685       18094       82102        3133           0

#Remove Extra columns
PDD_Target_Serovars_011025_Metadata <- select(PDD_Target_Serovars_011025_Metadata, -c(Extra1, Extra2))

#Upload list of all Enteritidis, Typhimurium, Infantis, Kentucky, Hadar, and Uganda isolates from NCBI ESearch used for Positive Homology Testing
Esearch_Target_Serovars_070323_Assemblies <- read.delim("all_target_assemblies_V2.txt", na.strings=c(""," ","NA"), col.names = "Assembly", header = FALSE)

#Remove the extra bits from Full.Assembly.Name 
Esearch_Target_Serovars_070323_Assemblies <- mutate(Esearch_Target_Serovars_070323_Assemblies,
                                                    Assembly = str_replace(Assembly, "_faas/Enteritidis_assemblies.txt:", ":"))

Esearch_Target_Serovars_070323_Assemblies <- mutate(Esearch_Target_Serovars_070323_Assemblies,
                                                    Assembly = str_replace(Assembly, "_faas/Hadar_assemblies.txt:", ":"))

Esearch_Target_Serovars_070323_Assemblies <- mutate(Esearch_Target_Serovars_070323_Assemblies,
                                                    Assembly = str_replace(Assembly, "_faas/Typhimurium_assemblies.txt:", ":"))

Esearch_Target_Serovars_070323_Assemblies <- mutate(Esearch_Target_Serovars_070323_Assemblies,
                                                    Assembly = str_replace(Assembly, "_faas/Uganda_assemblies.txt:", ":"))

Esearch_Target_Serovars_070323_Assemblies <- mutate(Esearch_Target_Serovars_070323_Assemblies,
                                                    Assembly = str_replace(Assembly, "_faas/Infantis_assemblies.txt:", ":"))

Esearch_Target_Serovars_070323_Assemblies <- mutate(Esearch_Target_Serovars_070323_Assemblies,
                                                    Assembly = str_replace(Assembly, "_faas/Kentucky_assemblies.txt:", ":"))

#Separate out Computed.types to create a IB_SeqSero2.serotype column
Esearch_Target_Serovars_070323_Assemblies <- separate(Esearch_Target_Serovars_070323_Assemblies, Assembly, sep = ":", into = c("Esearch_Serovar", "Assembly"), remove = FALSE)

#Get distribution
table(Esearch_Target_Serovars_070323_Assemblies$Esearch_Serovar, useNA = "always")
# Enteritidis       Hadar    Infantis    Kentucky Typhimurium      Uganda        <NA> 
#   37000        1704       14052       11100       25983         961           0 

#Copy to Excel
data.frame(table(Esearch_Target_Serovars_070323_Assemblies$Esearch_Serovar, useNA = "always")) %>%
  clipr::write_clip()


#Merge with PDD dataset
PDD_Esearch_Target_Serovars <- merge(PDD_Target_Serovars_011025_Metadata, Esearch_Target_Serovars_070323_Assemblies)
#83829

#Get distribution of Esearch assemblies with PDD information
table(PDD_Esearch_Target_Serovars$Esearch_Serovar, useNA = "always")
# Enteritidis       Hadar    Infantis    Kentucky Typhimurium      Uganda        <NA> 
#   34704        1636       13219       10800       22538         932           0 

#Copy to Excel
data.frame(table(PDD_Esearch_Target_Serovars$Esearch_Serovar, useNA = "always")) %>%
  clipr::write_clip()

#Determine distribution of isolates with SRA SRR Run IDs
table(is.na(PDD_Esearch_Target_Serovars$Run))
# FALSE  TRUE 
# 80073  3756 

#Merge with Pulsenet dataset
PDD_Esearch_PulseNet_Target_Serovars <- merge(PDD_Esearch_Target_Serovars, PulseNet_All_summ_CheckM_Pass, by.x = "Run", by.y = "SRR_id")
#1850

#Get distribution of Esearch assemblies with PDD information
table(PDD_Esearch_PulseNet_Target_Serovars$Esearch_Serovar, useNA = "always")
# Enteritidis       Hadar    Infantis    Kentucky Typhimurium      Uganda        <NA> 
#   896         181         361           3         382          27           0 