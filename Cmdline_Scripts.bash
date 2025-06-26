####NCBI Resources####

##Target Organism Assembly Website
#https://www.ncbi.nlm.nih.gov/assembly/GCF_000213635.1

##Target Organism Metadata
#Salmonella enterica subsp. enterica serovar Typhimurium str. UK-1 (enterobacteria)
#ASM21363v1
#Infraspecific name: Strain: UK-1
#BioSample: SAMN02602986
#BioProject: PRJNA63211
#Submitter: Arizona State University
#Date: 2011/05/13Assembly type: na
#Assembly level: Complete Genome
#Genome representation: full
#GenBank assembly accession: GCA_000213635.1 (latest)
#RefSeq assembly accession: GCF_000213635.1 (latest)
#RefSeq assembly and GenBank assembly identical: yes

##Manual downloading
#Go to website above and click Download Assembly
#Select either GenBank or RefSeq from Source database dropdown
#Select All file types (including assembly-structure directory) from File type dropdown
#Do this for both GenBank and RefSeq if you'd like, we did to allow for dual annotation later in process
#If you were to only choose one, do RefSeq since NCBI is adding GO Terms to their annotations

##Proteome used for this study
#GenBank Protein FASTA
#GCA_000213635.1_ASM21363v1_protein.faa

####Phase I - Subtractive Proteomics

###Vaxign2

##Running Vaxign2
#Run proteome of target through the Vaxign2 Dynamic Analysis
#https://violinet.org/vaxign2/

#Have several options to use such as FASTA Format or NCBI Protein, etc...do what works best for you

#The following options were used for our study:
#Select Pathogen Orgranism Type: Bacterium; Gram negative bacterium
#Checked Subcellular Localization (PSORTb)
#Checked Transmembrane Helix (TMHMM)
#Checked Adhesion Probability (SPAAN)
#Yes to Include Vaxign-ML Analysis
#No to Include Vaxitop Analysis

#Did batches of 250 since that seemed to be close to the limit of the server
#Once Vaxign is done running, click Copy to get all results and put into an txt file (Full_UK1_Vaxign2_Results.txt) in Phase_I_Subtractive_Proteomics folder
#Will result in 9 columns
#Have fun doing this a for all proteins

##Manipulating the Results for Easier Reading and Manipulation
#Copy the results over to an Excel file label tab "Vaxign2 Raw Results"
#Create a new tab, "Edited Vaxign2 Results" and copy over the raw results (make new tabs for each new step/filter below)
#Split the Localization(Probability) column into Localization and Localization Probability columns
#Did it in Excel via using "(" as the deliminator and removing extra bits from each columns (spaces, "Prob.=", and ")") via Ctrl + H
#Remove the extra spaces leftover in the Localization column as well
#I'm sure there is a way to do it via coding as well, have at it if you'd like :)

#Remove "Gene Accession", "Gene Symbol", "Locus Tag" columns since those did not populate
#Change Protein Accession to GenBank.Accession
#Change Protein Name to GenBank.Protein.Description
#Change rest of column names to being R friendly by adding "." or "_" instead of spaces

##Filtering Results
#You can do this via Excel pretty easily, but R code version is in Clean_Up.R (good way to just double check as well)
#If using R save the Edited Vaxign2 Results as Full_UK1_Vaxign2_Results.txt (SP_Output_1)
#Refer back to this file when prompted by the Clean_Up.R file

###Get Additional Annotations from Similar Salmonella Proteins

#The GenBank verion of UK1 is lacking in gene annotations, and Vaxign2 does not seem to recognize gene names
#Want to make a more informative table with original and additional annotations

#Copy just the Protein Accessions from the Excel and save in a separte tab-deliminated file (We used Notepad++)
#Name this file UK1_Ptn_Acc.txt (SP_Output_2)

#Log into SCINet and install and setup Miniconda

wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh

bash Miniconda3-py310_23.3.1-0-Linux-x86_64.sh

/project/fsepru113/dbradshaw/miniconda3

rm Miniconda3-py310_23.3.1-0-Linux-x86_64.sh

miniconda3/bin/conda update conda
miniconda3/bin/conda config --add channels defaults
miniconda3/bin/conda config --add channels bioconda
miniconda3/bin/conda config --add channels conda-forge

cd reverse_vaccinology/

mkdir UK1_ReAnno

cd UK1_ReAnno/

#Create an interactive session to run protein accessions from Vaxign2 in Blast
salloc -N1 -n4 -p short

#Get UK1 Accession Numbers to SCINet since blast remote does not work on local work computer
scp UK1_Ptn_Acc.txt david.bradshaw@ceres.scinet.usda.gov:/project/fsepru113/dbradshaw/reverse_vaccinology/UK1_ReAnno

#Activate blast conda environment
source /project/fsepru113/dbradshaw/miniconda3/bin/activate blast

#Run blast on accession numbers
#Limit blast to Salmonella both for project and to avoid Error: CPU usage limit was exceeded for remote blast
#https://www.biostars.org/p/9462196/#9487484

#Output the default output for future reference
blastp -query UK1_Ptn_Acc.txt -remote -db nr -out UK1_ReAnno_<date>.txt -evalue 1e-30 -max_target_seqs 10 -entrez_query "Salmonella [ORGN]"

#Output the archive output for future manipulation if needed
blastp -query UK1_Ptn_Acc.txt -remote -db nr -outfmt 11 -out UK1_ReAnno_archive_<date>.txt -evalue 1e-30 -max_target_seqs 10 -entrez_query "Salmonella [ORGN]"

#Output the targeted table from the archive file
blast_formatter -archive UK1_ReAnno_archive_<date>.txt -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -out UK1_ReAnno_tabular_<date>.txt

#Copy the information from SCINet down to you local computer into the appropriate folder by changing directories on your local computer to target location and doing scp
scp david.bradshaw@ceres.scinet.usda.gov:/project/fsepru113/dbradshaw/reverse_vaccinology/UK1_ReAnno/* .

#Outputs
#SP_Output_3 - UK1_ReAnno_<date>.txt
#SP_Output_4 - UK1_ReAnno_archive_<date>.txt
#SP_Output_5 - UK1_ReAnno_tabular_<date>.txt


#Open up UK1_ReAnno_tabular_<date>.txt in Excel, add the blast output names (qaccver, etc...) as column names, and replace any "N/A" in the stitle column with none [none] to help with filtering later
#Use Text to Columns on stitle column and deliminator [ to remove the Salmonella information, delete the subsequent taxonomic columns
#Make sure no extra bits are leftover i.e serovars with [ in name leading to extra columns
#Save as UK1_ReAnno_tabular_<date>_edited.txt (SP_Output_6)

#Upload into R to get gene names for the first subject Salmonella proteins with >98% identity or the first 500 hits to the query proteins

#Bye bye for now!

###Removal of Plasmid, Flagellar, and LPS Proteins
#Flagellar and LPS proteins are part of classification system for Salmonella serovars, likely to be too variable to be good epitope targets across serovars

#Create a new column called GenBank.Protein.Name
#Add the protein name from the GenBank.Protein.Description if it exists, otherwise add a " -"

#In Excel, create a new column called Category with following features
#If a protein has "plasmid" in its GenBank.Protein.Name or Similar.Salmonella.Proteins.Descriptions = Plasmid
#If a protein has "flagell*" in its GenBank.Protein.Name or Similar.Salmonella.Proteins.Descriptions = Flagellar
#If a protein has a named gene named, either GenBank or Similar Salmonella = Named
#If none of the above then = Unnamed

#Remove any proteins subsequently labled Plasmid or Flagellar

#Importing information from Excel...
###Non-flagellar, -LPS, and -plasmid UK1 Proteins = 167###

###Vaxijen Analysis

#Hope that went well for you, now onto antigenicity testing

#Use NCBI Genome Download to get fastas for all remaining proteins
#Save remaining GenBank protein accessions (n=167) as txt file (SP_Output_7) using nano
nano Pre_VaxiJen_UKI_GenBank_Proteins_Accessions.txt

#Run edirect targeting the protein database and saving the resulting fasta files
conda activate ncbi_edirect
cat Pre_VaxiJen_UKI_GenBank_Proteins_Accessions.txt | epost -db protein -format acc | efetch -format fasta > Pre_VaxiJen_UKI_GenBank_Proteins.fasta

#Run the resulting fasta file (SP_Output_8) in VaxiJen with threshold 0.5 and Summary Mode
#Copy and paste to Excel then sort based upon that column to remove spaces

#Removed extra information with Text to Column Fixed Width option if that seems ike it will work, beware of cutting stuff you need though

#For a more targeted approach:
#Find and Replace .1 with .1^, then Text to Columns to separate out Accession Numbers
#Text to Columns based upon ] to separate out Vaxijen Score and "(" to separate out description
#Find and Replace the following with blanks to remove extra bits: " )." ">", "[Salmonella enterica subsp. enterica serovar Typhimurium str. UK-1", and " Overall Protective Antigen Prediction = "

#Note on Excel, sometimes VaxiJen will have spaces left in number when using CTRL + H
#Copy and paste the space (which is a non-breaking white space) into CTRL + H to remove it
#https://superuser.com/questions/308072/how-to-remove-white-space-from-a-number

#Filter based upon a threshold of >=0.5
#https://pubmed.ncbi.nlm.nih.gov/34903179/
#https://pubmed.ncbi.nlm.nih.gov/36713875/

##Antigenic UK1 Proteins = 127

###Negative Homology Blasting

#Record information regarding number of proteins in NCBI for each potential host when you blasted them for future reference

##Gallus galus 
#Proteins - 070323 - 159,892 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9031/

##Meleagris gallopavo
#Proteins - 070323 - 30,663 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9103/

##Mus musculus
#Proteins - 070323 - 362,575 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/10090/

##Sus scrofa
#Proteins - 070323 - 98,748 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9823/

##Homo sapiens
#Proteins - 070323 - 1,907,751 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/

##Bos taurus
#Proteins - 070323 - 145,913 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9913/

##Blasting in SCINet

#Go to SCINet and add a folder to the reverse_vaccinology directory and go into it
mkdir UK1_Protein_Blasting
cd UK1_Protein_Blasting

#open up a new .txt file (SP_Output_9) and copy over GenBank Protein Accessions (n=127)
nano Post_VaxiJen_UK1_Protein_Accessions.txt

#In SCINet, open up a new .sh file with the target species name (SP_Output_Batch_1)
nano <species_name>_protein_negatve_blasting.sh

#In Notepad++, Open Draft_Negative_Homology_Protein_Blasting.sh (SP_Draft_Script_1) and use Ctrl + H to replace the species with the one you want to run
#Change the query in the blastp script as necessary as well
#Copy and paste over that information into your blank .sh file on SCINet

#Run the .sh file in a separate instances (SP_Batch_Output_2 & 3)
sbatch <species_name>_protein_negatve_blasting.sh



###Positive Homology Blasting

#Record information regarding number of proteins in NCBI for each serovar when you blasted them for future reference

#Salmonella - 070323 - 49,206,531 NCBI Proteins
#Salmonella enterica subsp. entrica - 43,234,599 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/59201/

#Infantis - 070323 - 963,345 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/595/

#Hadar - 070323 - 326,970 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/149385/

#Kentucky - 070323 - 1,123,476 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/192955/

#Uganda - 070323 - 94,853 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/487004/

#Enteritidis - 070323 - 4,241,000 NCBI Proteins
#https://www.ncbi.nlm.nih.gov/datasets/taxonomy/149539/

#In SCINet, open up a new .sh file with the target species name (SP_Output_Batch_4)
nano <species_name>_protein_negatve_blasting.sh

#In Notepad++, Open Draft_Positive_Homology_Protein_Blasting.sh (SP_Draft_Script_2) and use Ctrl + H to replace the species with the one you want to run
#Change the query in the blastp script as necessary as well
#Copy and paste over that information into your blank .sh file on SCINet

#Run the .sh file in a separate instances (SP_Output_Batch_5 & 6)
sbatch <species_name>_protein_negatve_blasting.sh

###Getting Protein Blasting information down to computer
#Copy the information from SCINet down to you local computer into the appropriate folder by changing directories on your local computer to target location and doing scp
scp david.bradshaw@ceres.scinet.usda.gov:/project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Protein_Blasting/*
#Delete what you believe you do not need (the stdout and stderr likely won't have much info for instance)

#Move over to R to upload all these files and then combine them with Post VaxiJen results in Excel

#...

#Importing information from Excel/R...
###Non-Similar to Hosts Proteins = 110###
###Similar to Other Salmonella Serovar Proteins = 101###

#...

###Adding Additional Annotations
#Save post-filtering proteins to a separate tab-deliminated file and upload into R to combine with other annotation information -> Post_Filtering_UK1_Proteins.txt

#See you back here for Phase II Epitope Discovery when you are done with Excel!

####PHASE II EPITOPE DISCOVERY

#Make a new directory 
mkdir Phase_II_Epitope_Discovery

#Copy and paste remaining protein accession numbers in a txt file using nano (ED_Output_1)
nano Post_Filtering_UK1_Protein_Accessions.txt 

#Download fasta for each of the remaining proteins (ED_Output_2)
epost -input Post_Filtering_UK1_Protein_Accessions.txt -db protein | efetch -format fasta > Post_Filtering_UK1.fasta

#Upload this fasta into NetMHCIpan for CTL/CD4/MHCI epitopes, NetMHCIIpan for HTL/CD8/MHCII epitopes, and BepiPred for LBL epitopes
#CTL - Cytotoxic T-Cell Lymphocyte
#HTL - Helper T-Cell Lymphocyte
#LBL - Linear B-Cell Lymphocyte

#Alleles used for this are based upon 

##NetMHCIpan
#https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
#v 4.1b
# Tmpdir made /var/www/html/services/NetMHCpan-4.1/tmp/netMHCpanblLW0J
# Input is in FSA format
# Peptide length 9
# Make both EL predictions
# Rank Threshold for Strong binding peptides   0.500
# Rank Threshold for Weak binding peptides   2.000
# Alleles: HLA-B40:06,HLA-B41:04,HLA-B41:03

#Copy and save the results to a new Excel for Phase II to manipulate it
#In the Excel, copy and paste header above to a file for future reference
#Use Text to Column Fixed width option on first column and split into appropriate columns
#Make sure it recognizes there is a BindLevel column
#Copy and paste to new sheet so you still have a full version, sort 2nd sheet by BindLevel by selecting columns, not just Ctrl + A
#Remove anything that is not a strong binder (SB)
#Adjust the Identity column to reflect actual protein name

#Unique_ID
#Create a Unique_ID column using the following formula
=CONCAT(K2,"_MHCI_",B2,"_Pos_",A2)
#Where the K column corresponds to Identity, B to MHC, and A to Pos

#...

#Importing information from Excel...
###885 Non-Unique MHCI Epitopes###

#...

##NetMHCIIpan 
#https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/
#v 4.0
#Input is FASTA format
#Peptide length is 15
#Prediciton Mode: EL
# Threshold for Strong binding peptides (%Rank) 0.5%
# Threshold for Weak binding peptides (%Rank)   2%
# Alleles: DRB1_1310,DRB1_1366,DRB1_1445,DRB1_1482

#Copy and save the results to the Phase II Excel to manipulate it
#In the Excel, copy and paste header above to a file for future reference
#Use Text to Column Fixed width option on first column and split into appropriate columns
#Make sure it recognizes there is a BindLevel column
#Copy and paste to new sheet so you still have a full version, sort 2nd sheet by BindLevel by selecting columns, not just Ctrl + A
#Remove anything that is not a strong binder (SB)

#Unique_ID
#Create a Unique_ID column using the following formula
=CONCAT(G2,"_MHCII_",B2,"_Pos_",A2)
#Where the G column corresponds to Identity, B to MHC, and A to Pos

#...

#Importing information from Excel...
###855 Non-Unique MHCII Epitopes###

#...

###Immunogenicity (MHCI ONLY)

##MHCI
#IEDB Analysis Resource 
#http://tools.iedb.org/immunogenicity/
#Just copy and paste the peptide sequences from the SB Excel Sheet, gotta love it when it is easy
#Copy results to MHCI SB Sheet, sort both tables by Peptide, and combine
#Remove those with a negative immunogenicity score

#...

#Importing information from Excel...
###396 Non-Unique Immunogenic MHCI Epitopes###

#...


###Antigencity (BOTH)

#Copy and paste Unique_ID and Peptide columns to txt files for MHCI and MHCII epitopes using nano

#MHCI
nano Pre_VaxiJen_UKI_MHCI_Epitope_Unique_Ids.txt
nano Pre_VaxiJen_UKI_MHCI_Epitope_Peptides.txt

#MHCII
nano Pre_VaxiJen_UKI_MHCII_Epitope_Unique_Ids.txt
nano Pre_VaxiJen_UKI_MHCII_Epitope_Peptides.txt

#Pop over to R really quick to make these into .fasta files

#...

#Howdy!

#Vaxijen
#http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html
#Copy results to MHCI Immunogenic Sheet, remove extra bits via Ctrl + H and Text to Column, sort both tables by Unique_Id, and combine
#Remove those with a VaxiJen score < 0.5

#...

#Importing information from Excel...
###210 Non-Unique MHCI Antigenic Epitopes###

#...

##MHCII
#Vaxijen
#http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html
#Copy results to MHCII SB Sheet, remove extra bits via Ctrl + H and Text to Column, sort both tables by Unique_Id, and combine
#Remove those with a VaxiJen score < 0.5

#...

#Importing information from Excel...
###446 Non-Unique MHCII Antigenic Epitopes###

#...

###Toxicity and Hydropathicity
#ToxinPred
#https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.php
#Copy and paste peptide sequences into ToxinPred, let it run with defaults, and copy down results

#First make a sheet that just removes the toxic epitopes for each MHC Class

#...

#Importing information from Excel...
###210 Non-Unique MHCI Non-Toxic Epitopes###
###440 Non-Unique MHCII NOn-Toxic Epitopes###

#...

#Then make a sheet that further removes the hydrophobic epitopes from each MHC Class
#Hydropathicity corresponds with the GRAVY Index from ProtParam
#A negative value is a hydrophilic sequence
#Remove those with a positive value

#...

#Importing information from Excel...
###132 Non-Unique MHCI Hydrophilic Epitopes###
###394 Non-Unique MHCII Hydrophilic Epitopes###

#...

###Epitope Genome Blasting Per Serovar

#Since this blasting will be a combination of MHCI and MHCII epitopes and we are going to need this data anyways for BepiPred, save the Pre Positive Homology Results as tab-deliminated txts and upload into R
#Outputs
#ED_Output_9 - Pre_Positive_Homology_UK1_MHCII_Results.txt
#ED_Output_10 - Pre_Positive_Homology_UK1_MHCI_Results.txt
#ED_Output_11 - Pre_Positive_Homology_UKI_MHCI_Epitopes.fasta
#ED_Output_12 - Pre_Positive_Homology_UKI_MHCII_Epitopes.fasta

#...


#In a local terminal, combine the two fasta into one file (ED_Output_13_
cat Pre_Positive_Homology_UKI_MHCI_Epitopes.fasta Pre_Positive_Homology_UKI_MHCII_Epitopes.fasta > Pre_Positive_Homology_UKI_Epitopes.fasta

#Go to SCINet and add a folder to the reverse_vaccinology directory and go into it
mkdir UK1_Epitope_Proteome_Blasting
cd UK1_Epitope_Proteome_Blasting

#Get the combined epitope fasta to SCINet
scp Pre_Positive_Homology_UKI_Epitopes.fasta david.bradshaw@ceres.scinet.usda.gov:/project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting

#Do the following steps for each of the seorvars being tested
#Salmonella enterica subsp. enterica serovar Infantis (taxid:595)
#Salmonella enterica subsp. enterica serovar Hadar (taxid:149385)
#Salmonella enterica subsp. enterica serovar Uganda (taxid:487004)
#Salmonella enterica subsp. enterica serovar Enteritidis (taxid:149539)
#Salmonella enterica subsp. enterica serovar Kentucky (taxid:192955)
#Salmonella enterica subsp. enterica serovar Typhimurium (taxid:90371)

#Make a directory for each of the serovars
mkdir <serovar>_faas

#Move/Copy the Pre_Positive_Homology_UK1_Epitopes.fasta file to each serovar's working directory
cp Pre_Positive_Homology_UKI_Epitopes.fasta <serovar>_faas

#After cding into into each serovar's directory, make a list of all of the possible urls that could be downloaded (ED_Output_Batch_1)
esearch -db assembly -query txid192955[organism:exp] | esummary   | xtract -pattern DocumentSummary -element FtpPath_GenBank   | while read -r url ; do      path=$(echo $url | perl -pe 's/(GC[FA]_\d+.*)/\1\/\1_protein.faa.gz/g') ; echo $path ; done > Kentucky_urls.txt

#Determine and record what that number is
wc <serovar>_urls.txt -l

#Enteritidis - 37433
#Hadar - 1710
#Infantis - 14079
#Kentucky - 11118
#Typhimurium - 28136
#Uganda - 961

#Actually attempt to download all of the faa protein files from NCBI by editing the Proteome_Downloading_Draft.sh in text editor and copy to nanoed version in SCINet
nano <serovar>_proteomes.sh
sbatch <serovar>_proteomes.sh

#Once finished, determine and record the number of initially downloadable faa files
ls protein_data | wc -l
#Enteritidis - 36831
#Hadar - 1685
#Infantis - 13974
#Kentucky - 11016
#Typhimurium - 25549
#Uganda - 949

#Decompress all of the files
gzip -d *.gz

#Gzip will fail if it comes across a corrupted file, thus mv it out of the protein_data directory into the <serovar>_faas directory and restart decompression
mv <corrupted_file> ../
gzip -d *.gz

#Once you are done with initial set of decompression, determine and record the number of correctly downloaded faa files
ls protein_data | wc -l
#Enteritidis - 36804
#Hadar - 1685
#Infantis - 13953
#Kentucky - 11016
#Typhimurium - 25029
#Uganda - 949

##Dealing with Corrupted Files

#Try redownloading the corrupted files by first, moving back a directory
cd ../

#Determine and record the number of corrupted files
ls *.gz | wc -l
#Enteritidis - 27
#Hadar - 0
#Infantis - 21
#Kentucky - 0
#Typhimurium - 520* (Assumed rest was corrupted when gzip kept hitting corrupted files)
#Uganda - 0

#Make a new directory to act as an intermediatary
mkdir missing_faas

#Use a for loop to take the names of the corrupted files, search the urls, and try to redownload the files
for i in ls *.gz ; do path=$(grep "$i" <serovar>_urls.txt) ; wget echo $path -P missing_faas ; done

#Remove the extra htmls
rm missing_faas/index*

#Attempt to decompress these files
gzip -d missing_faas/*

#If sucessful, move them to the main proteome location
mv missing_faas/* protein_data/

#Remove the corrupted files
rm *.gz

##Retry GenBank Downloading

#Make a list of all of the downladed assemblies
ls protein_data/ | cut -d "_" -f 1,2 > <serovar>_downloaded_GB_proteomes_Rd1.txt

#Get the assemblies from all the possible downloads
cut <serovar>_urls.txt -d "/" -f 10 | cut -d "_" -f 1,2 > <serovar>_ftp_assemblies.txt

#Find the assemblies that are in the possible urls but were not downloaded
grep -Fxv -f <serovar>_downloaded_GB_proteomes_Rd1.txt <serovar>_ftp_assemblies.txt > <serovar>_missing_proteomes_post_GB_Rd1.txt

#Determine and record the number of missing assemblies
wc <serovar>_missing_proteomes_post_GB_Rd1.txt -l
#Enteritidis - 575
#Hadar - 25
#Infantis - 105
#Kentucky - 102
#Typhimurium - 2587
#Uganda - 12

#Attempt to download the GenBank versions again
for i in $(cat <serovar>_missing_proteomes_post_GB_Rd1.txt); do GenBankpath=$(grep "$i" <serovar>_urls.txt); wget -q --show-progress "$GenBankpath" -P missing_faas;  done

#If sucessful, determine and record the number of missing GenBank assemblies
ls missing_faas | wc -l
#Enteritidis - 0
#Hadar - 1
#Infantis - 10
#Kentucky - 2
#Typhimurium - 9
#Uganda - 1

#If sucessful, attempt to decompress them
gzip -d missing_faas/*.gz

#If sucessful, move them to the main proteome location
mv missing_faas/* protein_data/

##RefSeq Assembly Download

#Make a new downloaded assemblies file
ls protein_data/ | cut -d "_" -f 1,2 > <serovar>_downloaded_GB_proteomes_Rd2.txt

#Make a new missing assemblies file
grep -Fxv -f <serovar>_downloaded_GB_proteomes_Rd2.txt <serovar>_ftp_assemblies.txt > <serovar>_missing_proteomes_post_GB_Rd2.txt

#Attempt to download the RefSeq version of these assemblies
for i in $(cat <serovar>_missing_proteomes_post_GB_Rd2.txt); do GenBankpath=$(grep "$i" <serovar>_urls.txt); RefSeqpath=$(echo $GenBankpath | sed -r 's/GCA/GCF/g') ; wget -q --show-progress "$RefSeqpath" -P missing_faas;  done

#If sucessful, determine and record the number of RefSeq GCF assemblies
ls missing_faas | wc -l
#Enteritidis - 142
#Hadar - 18
#Infantis - 68
#Kentucky - 82
#Typhimurium - 425
#Uganda - 11

#If sucessful, attempt to decompress them
gzip -d missing_faas/*.gz

#If sucessful, move them to the main proteome location
mv missing_faas/* protein_data/

#Determine and record the number of proteomes
ls protein_data/ | wc -l
#Enteritidis - 37000
#Hadar - 1704
#Infantis - 14052
#Kentucky - 11100
#Typhimurium - 25983
#Uganda - 961

#Make a final downloaded assemblies file
ls protein_data/ | cut -d "_" -f 1,2 > <serovar>_downloaded_GB_RF_proteomes.txt

#For the purposes of finding missing assemblies, convert all GCF to GCA
sed  -i 's/GCF/GCA/g' <serovar>_downloaded_GB_RF_proteomes.txt

#Make a final missing assemblies file
grep -Fxv -f <serovar>_downloaded_GB_RF_proteomes.txt <serovar>_ftp_assemblies.txt > <serovar>_missing_proteomes_final.txt

#Double check a few of these files to make sure they do not have .faa files by going to the NCBI Assemblies search page
https://www.ncbi.nlm.nih.gov/assembly

#Determine and record the number of missing proteomes
wc <serovar>_missing_proteomes_final.txt -l
#Enteritidis - 433
#Hadar - 6
#Infantis - 27
#Kentucky - 18
#Typhimurium - 2153
#Uganda - 0

#Remake a final downloaded assemblies file
ls protein_data/ | cut -d "_" -f 1,2 > <serovar>_downloaded_GB_RF_proteomes.txt

#Remove the temporary missing_faas folder
rm -rf missing_faas/

#Make directories to hold results in each serovar's working directory
mkdir blast_dbs
mkdir epitope_blasting_results
mkdir slurm_stds

#Adjust the Draft_Positive_Homology_Epitope_Proteomes_Blasting.sh in Notepad++ to your serovar of interest and run a test via sbatch array
#Set Line 5 to 1-2 to test out your set up
nano test.sh
sbatch test.sh

#Double check the output folders and files for correct info
#If everythign is good, remove the test files
rm blast_dbs/*
rm epitope_blasting_results/*
rm std*

#Adjust the array to 1 + the number of faa files if below 10000 (x)
#Line 5 -> SBATCH --array=1-(1+x)
#etc...
ls protein_data/ | wc -l
nano <serovar>_proteomes_epitope_blasting.sh

#If above 10000, then use the Draft_Positive_Homology_Epitope_Proteomes_Blasting.sh for first 10000
#Then adjust the Draft_Positive_Homology_Epitope_Proteomes_Blasting_Pt2.sh for the rest
#Keep in mind to adjust the count accordingly in line 14 for each new 10000
#Line 5 -> #SBATCH --array=1-4062
#Line 14 -> proteome_file="$(ls protein_data/*.faa | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#Adjustments above are acually running proteomes 10001-14062

#Let the array do its thing
sbatch <serovar>_proteomes_epitope_blasting.sh

#Move all std* files to slurm_stds directory
mv std* slurm_std/

##One-step version of summarizing blast results(OSSB)  (slowers, less hands on, possible to use for <10k proteomes)

#Combine all blast results into one file
cd epitope_blasting_results
grep "" *.txt > <serovar>_proteomes_epitope_blasting_results.tsv
mv <serovar>_proteomes_epitope_blasting_results.tsv ../<serovar>_proteomes_epitope_blasting_results_<date_of_url_list_download>.txt

#Make a slimmed down version of this for R
cat Uganda_proteomes_epitope_blasting_results_070723.txt  | cut -f 1,3  | sed 's/_tabular_070723.txt:/;/g' > Uganda_proteomes_epitope_blasting_results_070723_R.txt

#Combine all stdout and stderr into one file
cd ../
cat slurm_std/stderr* > <serovar>_proteomes_epitope_blasting_stderr.txt
cat slurm_std/stdout* > <serovar>_proteomes_epitope_blasting_stdout.txt

#The resulting <serovar>_proteomes_epitope_blasting_results_<date_of_url_list_download>_R.txt will likely be very large thus want to use lots of cores in SCINet to run in R

# Adjust the Draft_Proteome_Blasting_Summaries_Parallel_Processing_R.sh for your serovar of interest and for the cores you wish to run and copy over to a nanoed file
nano <serovar>_blast_results_summarizing.R

#Adjust the Draft_Proteome_Blasting_Summaries_Parallel_Processing.sh for your serovar fo interest and for the cores (-c) you wish to run and copy over to a nanoed file
#Adjust the --mem-per-cpu and -p (node type) to the size and time you want
#Example of setups that worked
#Kentucky - JobId 9996182; Input file size = 1.8G; -c 12; --mem-per-cpu 12gb; -p mem; time took to run = 4-17:27:20
#Infantis - JobId 10005519; Input file size = 2.3G; -c 36; --mem-per-cpu 12gb; -p mem; time took to run = 4-08:08:59
nano <serovar>_blast_results_summarizing.sh
sbatch <serovar>_blast_results_summarizing.sh


#Copy the following files from SCINet down to you local computer into the appropriate folder by changing directories on your local computer to target location and doing scp
#<serovar>_proteomes_epitope_homology_counts_summary.txt
#<serovar>_proteomes_epitope_homology_summary.txt
#<serovar>_proteomes_epitope_homology.txt
#<serovar>.RData

#Delete what you believe you do not need (the stdout and stderr likely won't have much info for instance)

#Head over to R for summarization of the number of epitopes with 100% identity with varying percentages of proteomes per serovar

##Multi-step array version of summarizing blast results (MSASB) (quicker, more hands on, definitely use for >10k proteomes)

#Make a directory to hold per proteome summaries from the array
mkdir <serovar>_faas/epitope_blasting_summaries_Step_1/

#Edit the Draft_Proteome_Blasting_Summaries_Step_1_R.sh to reflect serovar
nano <serovar>_proteome_blasting_summaries_Step_1.R

#Edit the Draft_Proteome_Blasting_Summaries_Step_1_Pt1.sh to reflect serovar
#Set Line 4 to 1-2 to test out your set up
#Copy over to nanoed file and let it run
nano Enteritidis_proteome_blasting_summaries_Pt1.sh
sbatch Enteritidis_proteome_blasting_summaries_Pt1.sh

#Double check the output folders and files for correct info
#If everythign is good, remove the test files
rm <serovar>_faas/epitope_blasting_summaries
rm std*

#Adjust the array to the number of faa files if below 10000 (x) and run it
#Line 4 -> SBATCH --array=1-x
ls protein_data/ | wc -l
nano Enteritidis_proteome_blasting_summaries_Pt1.sh
sbatch Enteritidis_proteome_blasting_summaries_Pt1.sh

#If above 10000, then use the Draft_Proteome_Blasting_Summaries_Step_1_Pt1.sh for first 10000
#Then adjust the Draft_Proteome_Blasting_Summaries_Step_1_Pt2-x.sh for the rest
#Keep in mind to adjust the count accordingly in line 12 for each new 10000
#Line 5 -> #SBATCH --array=1-4062
#Line 14 -> proteome_file="$(ls protein_data/*.faa | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#Adjustments above are acually running proteomes 10001-14062
nano Enteritidis_proteome_blasting_summaries_Ptx.sh

#Let the array do its thing
sbatch Enteritidis_proteome_blasting_summaries_Ptx.sh

#Make a new directory and move all the slurm stderrs and stdouts to it
mkdir Enteritidis_faas/epitope_blasting_summaries_Step_1
mv stderr* Enteritidis_faas/epitope_blasting_summaries_Step_1/
mv stdout* Enteritidis_faas/epitope_blasting_summaries_Step_1/


#Determine what proteomes were actually summarized and record the number of proteomes
ls <serovar>_faas/epitope_blasting_summaries/ | wc -l
#Enteritidis - 26985
#Typhimurium - 25957 

#Make a list of the summarized proteomes
ls <serovar>_faas/epitope_blasting_summaries/ | cut -d "_" -f 1,2 > <serovar>_summaries_Rd1.txt

#Make a missing summarized proteomes file and record the number
grep -Fxv -f <serovar>_summaries.txt <serovar>_downloaded_GB_RF_proteomes.txt > <serovar>_missing_proteomes_post_summaries_Rd1.txt
#Enteritidis - 10015
#Typhimurium - 26

#Make a list of blast results files
ls <serovar>_faas/epitope_blasting_results > <serovar>_blasting_results_files.txt
#Enteritidis - 37000
#Typhimurium - 25983

#Make a temporary directory
mkdir <serovar>_missing_summaries

#Copy blasting results that were not summarized to this folder
for i in $(cat Typhimurium_missing_proteomes_post_summarizing_Rd1.txt); do Results=$(grep "$i" Typhimurium_blasting_results_files.txt); Results_path="Typhimurium_faas/epitope_blasting_results/$Results"; cp $Results_path Typhimurium_missing_summaries/; done

#Edit the Draft_Proteome_Blasting_Summaries_Step_1_Redo_Pt1.sh to reflect serovar and number of missing summaries
#Copy over to nanoed file and let it run
nano Enteritidis_proteome_blasting_summaries_Step_1_redo_Pt1.sh
sbatch Enteritidis_proteome_blasting_summaries_Step_1_redo_Pt1.sh

#If need to do more than 10000, then use same strategy as above and adjust Draft_Proteome_Blasting_Summaries_Step_1_Redo_Pt2-x.sh

#Double check you now have all the necessary summaries and record number
ls Typhimurium_faas/epitope_blasting_summaries/ | wc -l
#Enteritidis - 37000
#Typhimurium - 25983

#If you do then cd into that directory, if not then redo above steps
cd Typhimurium_faas/epitope_blasting_summaries/

#Get the name of the first file in the folder and use in next script
ls * | head -n 1

#Concatenate all summaries by first grabbing header from one and then results from all
head -n 1 GCA_000006945.1_proteome_epitope_homology.txt > Typhimurium_proteomes_epitope_homology.tsv; tail -n +2 -q *.txt >> Typhimurium_proteomes_epitope_homology.tsv

#If need be, you can also use edit and use Draft_Proteome_Blasting_Summaries_Post_Step_1_Concatenation.sh

#Number of lines should be about number of epitopes x number of downloaded proteomes
#Theoretical Enteritidis - 526 x 37000 = 19462000
#Theoretical Hadar - 526 x 1704 = 896304
#Theoretical Infantis - 526 x 14052 = 7391352
#Theoretical Kentucky - 526 x 11100 = 5838600
#Theoretocal Typhimurium - 526 x 25983 = 13667058
#Theoretical Uganda - 526 x 961 = 505486

wc -l Typhimurium_proteomes_epitope_homology.tsv

#Actual Enteritidis - 19461855
#Actual Hadar - 896304
#Actual Infantis - 7391322
#Actual Kentucky - 5838552
#Actual Typhimurium - 13667053
#Actual Uganda - 505486

#Move the remaining stdouts and stderrs to previous folder
mv stderr* Enteritidis_faas/epitope_blasting_summaries_Step_1/
mv stdout* Enteritidis_faas/epitope_blasting_summaries_Step_1/

#Move to the main directory for next steps in R
mv Typhimurium_proteomes_epitope_homology.tsv ../../<serovar>_proteomes_epitope_homology.txt

#Make a new directory to hold per epitope splitting of <serovar>_proteomes_epitope_homology.txt
mkdir <serovar>_splitting

Change into this directory and copy over <serovar>_proteomes_epitope_homology.txt
cd <serovar>_splitting
cp ../<serovar>_proteomes_epitope_homology.txt .

#Split <serovar>_proteomes_epitope_homology.txt into a summary of each epitope (Unique_Id column = column 1)
awk '{print > $1".txt"}' <serovar>_proteomes_epitope_homology.txt

#Remove the extra files
rm <serovar>_proteomes_epitope_homology.txt
rm Unique_Id.txt

#Check to make sure there is the correct number of files i.e. matches number of non-unique epitopes
ls * | wc -l #526

#Make a directory to hold results of Draft_Proteome_Blasting_Summaries_Step_2
mkdir Infantis_faas/epitope_blasting_summaries_Step_2

Enteritidis_proteome_blasting_summaries_Part_2.R

#Adjust Draft_Proteome_Blasting_Summaries_Step_2.sh draft_blast_results_summarizing_Part_2_array.sh to your serovar and nanofy it
#start wtih 1-2 to test setup
nano <serovar>_Enteritidis_proteome_blasting_summaries_Part_2.sh
sbatch <serovar>_Enteritidis_proteome_blasting_summaries_Part_2.sh 

#Check that the output is expected and remove test files
cat Infantis_faas/epitope_blasting_summaries_Step_2/*
rm Infantis_faas/epitope_blasting_summaries_Step_2/*
rm std*

#Adjust the sh file to the total number of epitopes tested (n=526) and let it run
nano <serovar>_Enteritidis_proteome_blasting_summaries_Part_2.sh
sbatch <serovar>_Enteritidis_proteome_blasting_summaries_Part_2.sh 

#Move into the directory and combine all files into one
cd Infantis_faas/epitope_blasting_summaries_Step_2

head -n 1 AEF05957.1_MHCI_HLA-B*40:06_Pos_35_proteomes_epitope_homology_summary.txt > Infantis_proteomes_epitope_homology_summary.tsv; tail -n +2 -q *.txt >> Infantis_proteomes_epitope_homology_summary.tsv

#Move it to the main directory
mv Infantis_proteomes_epitope_homology_summary.tsv ../../Infantis_proteomes_epitope_homology_summary.txt

#Make a directory to store the stderr and stdout files and move them there
mkdir Infantis_faas/epitope_blasting_summaries_Step_2_stds
mv std* Infantis_faas/epitope_blasting_summaries_Step_2_stds


#Copy the following files from SCINet down to you local computer into the appropriate folder by changing directories on your local computer to target location and doing scp
#<serovar>_proteomes_epitope_homology_summary.txt
#<serovar>_proteomes_epitope_homology.txt
#Delete what you believe you do not need (the stdout and stderr likely won't have much info for instance)

#Head over to R for summarization of the number of epitopes with 100% identity with varying percentages of proteomes per serovar

#...

#Importing information from R & Excel...
###69 Non-Unique MHCI with 100% Identity to 99% of All Serovars' Proteomes###
###216 Non-Unique MHCII with 100% Identity to 99% of All Serovars' Proteomes###

#...

###BepiPred Linear B-Cell Epitope Identification
#https://services.healthtech.dtu.dk/services/BepiPred-3.0/

#Can only run 50 proteins at a time so have to run in batches of 50 from the fasta
#We ran all the proteins so that we have this information on hand if need be
#NOTE: If protein is larger than what BepiPred can handle, split it in half and label "Part 1" and "Part 2" after protein description
#Example: AEF09523.1

##Check the following options: 
#For "Top epitope percentage cutoff", select "Higher confidence (top 20%)"
#Run with default "Threshold for predicting B-cell epitope residues" = 0.1512
#Select "Yes" for "Use sequential smoothing (lineaer epitope prediction mode) on B-cell epitope probability score graphs (see instructions):"

##Post-Processing:
#Click "Download bebipred3 results (zip)"
#Extract out the raw_output
#Combine all results onto one Excel sheet in epitope summary excel
#Remove the extra bits in the Accession column by replacing .1 with .1^ and then using Text to Column to separate by ^
#Relabel columns to the following: Accession, Residue, Score, Rolling.Mean.Score
#Save this sheet as its own tab-deliminated txt file -> BepiPred_Raw_Results.txt

#Head over to R for lots of loops and branching filtering of epitopes and then finding unqiue epitopes

###BepiPred Linear B-Cell Epitope Identification

#...

###Find Unique Epitopes

#...

#Importing information from R...
###42 Unique MHCI with 100% Identity to 99% of All Serovars' Proteomes###
###126 Unique MHCII with 100% Identity to 99% of All Serovars' Proteomes###

#...

###Find Nonoverlapping, Unique Epitopes

#Once you have copied down the Post_Filtering_UK1_Peptide_Summary dataframe, you need to manipulate it to find and combine overlapping epitopes

#Sort by Identity (GenBank Accession) and and then Pos (descending) 
#Create a temporary column called Order and poplate with numbers designating the order of the rows to help with analysis
#Visually find any epitopes that seem like they could be overlapping per their position within the protein
#It may be helpful to delineate proteins by using thicker borders 
#Once you find overlapping sections, determine which one has the highest VaxiJen_Score and bold it
#It may be helpful to color regions of overlapping with alternating colors to keep track of things

#Copy the full table over to new sheet and delete any colored row that is not bolded


#Create three new colums after the Bcell_epitope_percentage column: LBL_Epitope, GT2_Alleles, GT1_Allele

#LBL_Epitope
#Sort by MHC_Type and then the Bcell_epitope_percentage to determine what T-Cell epitopes are also B-Cell epitopes
#If it is MHCI, then want Bcell_epitope_percentage >= 88% AND 8 or more continous LBL residues
#If it is MHCII, then want Bcell_epitope_percentage >=53.3%  AND 8 or more continous LBL residues
#If it is a LBL, then label it TRUE, otherwise FALSE

#GT2_Alleles or GT1_Allele
#Look at the MHC_Alleles column and if it has >=2 Alleles label it as TRUE for GT1_Allele, otherwise FALSE
#Look at the MHC_Alleles column and if it has >=3 Alleles label it as TRUE for GT2_Allele, otherwise FALSE


#Sort by Order, and delete the Order column from both sheets (Unique and NonOverlapping)

#Save sheet as Post_Filtering_UK1_Nonoverlapping_Unique_Peptide_Summary.txt

#...

#Importing information from Excel...
###41 Unique, Nonoverlapping MHCI with 100% Identity to 99% of All Serovars' Proteomes###
###52 Unique, Nonoverlapping MHCII with 100% Identity to 99% of All Serovars' Proteomes###

#...

###Summarizing Proteins Associated with  Final Epitope Table

#Go to R to upload the nonoverlapping, unique table and make an protein summary and annotated version of it

#...
#Importing information from R and Excel...
###34 proteins with MHCI epitopes with 100% Identity to 99% of All Serovars' Proteomes###
###39 proteins with MHCII epitopes with 100% Identity to 99% of All Serovars' Proteomes###
###53 proteins with either MHCI and/or MHCII epitopes with 100% Identity to 99% of All Serovars' Proteomes###
###20 proteins with both MHCI and MHCII epitopes with 100% Identity to 99% of All Serovars' Proteomes###

#...

####Phase III Construct Construction

#Use the R script to create constructs (mulltiepitope, peptide, polypeptide) with the remaining epitopes 

#...

###Multiepitope Construct - Most Antigenic Epitopes (9) Localization Agnostic x MHC Type Design

#Importing information from R and Excel...

###Poultry Linkers Epitope Region - 441 aa ###
###EAAAKIEGEDMRLAAAYGEDRRTLNVAAYTEREGKAAAAAYKEDNELREAAAYSEADVQGHVAAYYEYNFRTAYAAYYEKTDNTRMAAYKDKAFDVKLAAYKETGERLSIGPGPGASGDLTVEVKESDGSGPGPGAQKLAIEIRDGDQRRGPGPGRAGYRADVKNNDSNVGPGPGLHYFSDDKGSDGDQTGPGPGQNIAVVRRADGSGTSGPGPGHWEITNTFRYRINEHGPGPGTPGLRFDHHSIVGDNGPGPGQGNPVTGTDKQAVSPGPGPGTALTFSRDGKTQDKNKKLENEFKGRAKKTWHARFAYDKEKTDRKKKRPFAGNTGTVDDKDKKAAYSNSKRTNDQQDRKKSRGNYRYTDKDLVKYKKARYRFEYVRRSSDIRKKSDGTKINYANKVINNKKDERVALREAKKAENATTDKAKKKERIAEKGAEAAAK###

#...

####Phase IV Construct Evaluation

#Copy Flagelin and Epitope Region sequences to text editor to save a construct fasta file

##Constructs tested
#Any_Loc_Multiepitope_Chicken_Linkers_GFP_Detached

#Create an Excel Sheet to store all information

#For all the tools below, add a Job Name and Email if offered

###General Construct Properties
#Copy all results from this section to the Excel for summarization

##Antigencity

#VaxiJen v.2.0
#Upload the construct fasta file to server
#Target Organism is Bacteria
#Check Summary mode
#Threshold = 0.5
#http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html

#ANTIGENpro
#SCRATCH Protein Predictor has option to just do ANTIGENpro as same time as SOLpro so just did it as another check
#Copy protein sequence (no > line) separately
#http://scratch.proteomics.ics.uci.edu/

##Allergenicity

#AllergenFP v.1.0
#Copy protein sequence (no > line) separately
#http://ddg-pharmfac.net/AllergenFP/index.html

#AllerTOP v.2.0
#Copy protein sequence (no > line) separately
#https://www.ddg-pharmfac.net/AllerTOP/index.html

##Toxicity
#ToxinPred2
#Upload the construct fasta file to server
#Machine Learning Technique used for developing model - Hybrid (RF+BLAST+MERCI) is Default
#https://webs.iiitd.edu.in/raghava/toxinpred2/batch.html

##Physiochemial Properties
#Expasy-ProtParam
#Copy protein sequence (no > line) separately
#https://web.expasy.org/protparam/

##Solubility
#SOLpro
#Copy protein sequence (no > line) separately
#http://scratch.proteomics.ics.uci.edu/

###Structure Assessment

##Secondary Structure Prediction
#PSIPRED v4.0
#Copy protein sequence (no > line) separately
#http://bioinf.cs.ucl.ac.uk/psipred/

#Summarizing predictions
#Save the PNG of the Sequence Plot for Supplementary/Manuscript Figures
#Save the SS2 Format Output and open in Excel
#Use Unique on Column C to get unique occurances of seconday structure types
#Use Countif in cell next to each unique secondary structure type to test the same range for the occurances of that secondary structure Type
#Next to the Countif column make the percentage Column
#Save the SS2 format as an Excel file
#Copy and paste summary to Phase IV Excel

##Tertiary Structure Prediction
#Phyre2
#Copy protein sequence (no > line) separately, batch processing only allowed on Normal mode
#Run in Intensive mode
#http://www.sbg.bio.ic.ac.uk/phyre2/html/page.cgi?id=index

##Tertiary Structure Refinement
#GalaxyRefine
#Although GalaxyRefine2 exists, it has a 300 AA limit
#Copy protein sequence (no > line) separately
#https://galaxy.seoklab.org/cgi-bin/submit.cgi?type=REFINE

#Choose the best refined model using the following method
#Copy the table down to Excel
#Assign 1-5 (1 being best and 5 being worst) for the three main parameters in the following manner:
#GDT-HA = 1 = highest, 5 = lowest
#RMSD = 1 = lowest, 5 = highest
#MolProbity = 1 = lowest, 5 = highest
#Ties are given the same better number
#Find the sum of all three values for all five models, and choose the lowest sum
#Ties go to model with lowest MolProbity Score

##Tertiary Structure Validation

#PDBSum
#Upload original Phyre2 and GalaxyRefine chosen refined model separately
#Download the PROCHECK Figure and print to save as pdf the results pages

#https://saves.mbi.ucla.edu/
#ProSA-Web
#Upload original Phyre2 and GalaxyRefine chosen refined model separately

##CBL prediction
#ElliPro
#Upload chosen refined model 
##Molecular Docking


##Molecular Docking Residue Interactions Validation

#Search for your reference files in the PDBsum section (e.g. 3v47 - zebrafish TLR5 homodimer with flagellin)
#Once there click on the Interface Summary or Summaries you are interested on the left side (e.g. A-C, B-C, A-D, B-D for 3v47)
#Then scroll down to the "Residue interactions across interface" and open the "List of interactions" in a new tab
#Copy and paste all of this to a new txt file and save it with a unique Name

#3v47-Flagellins
#Create total interactions residues list for each chain in the reference
grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" 3v47_TLR5_Complex_Chains_A_C_Interactions.txt | while read p; do echo $p; done | cut -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > 3v47_TLR5_Complex_Chains_A_C_Interactions_A_Residues.txt
grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" 3v47_TLR5_Complex_Chains_A_D_Interactions.txt  | while read p; do echo $p; done | cut -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > 3v47_TLR5_Complex_Chains_A_D_Interactions_A_Residues.txt
grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" 3v47_TLR5_Complex_Chains_B_C_Interactions.txt  | while read p; do echo $p; done | cut -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > 3v47_TLR5_Complex_Chains_B_C_Interactions_A_Residues.txt
grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" 3v47_TLR5_Complex_Chains_B_D_Interactions.txt  | while read p; do echo $p; done | cut -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > 3v47_TLR5_Complex_Chains_B_D_Interactions_A_Residues.txt

#Create total interactions residues list for each chain in the reference
grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" 3v47_TLR5_Complex_Chains_A_C_Interactions.txt  | while read p; do echo $p; done > 3v47_TLR5_Complex_Chains_A_C_Interactions_Only_SS.txt
mapfile -t my_array < <( grep -n "^1. " 3v47_TLR5_Complex_Chains_A_C_Interactions_Only_SS.txt | cut -f1 -d: )
if [ ${#my_array[@]} == 3 ]; then rm1="${my_array[1]}"; rm2=$(echo "${my_array[2]}"-1 |bc); sed -e "$rm1,$rm2"'d' 3v47_TLR5_Complex_Chains_A_C_Interactions_Only_SS.txt > 3v47_TLR5_Complex_Chains_A_C_Bonds_Only_SS.txt; else rm1="${my_array[1]}"; sed -e "$rm1,$"'d' 3v47_TLR5_Complex_Chains_A_C_Interactions_Only_SS.txt > 3v47_TLR5_Complex_Chains_A_C_Bonds_Only_SS.txt; fi
cut 3v47_TLR5_Complex_Chains_A_C_Bonds_Only_SS.txt -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > 3v47_TLR5_Complex_Chains_A_C_Bonds_A_Residues.txt

grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" 3v47_TLR5_Complex_Chains_A_D_Interactions.txt  | while read p; do echo $p; done > 3v47_TLR5_Complex_Chains_A_D_Interactions_Only_SS.txt
mapfile -t my_array < <( grep -n "^1. " 3v47_TLR5_Complex_Chains_A_D_Interactions_Only_SS.txt | cut -f1 -d: )
if [ ${#my_array[@]} == 3 ]; then rm1="${my_array[1]}"; rm2=$(echo "${my_array[2]}"-1 |bc); sed -e "$rm1,$rm2"'d' 3v47_TLR5_Complex_Chains_A_D_Interactions_Only_SS.txt > 3v47_TLR5_Complex_Chains_A_D_Bonds_Only_SS.txt; else rm1="${my_array[1]}"; sed -e "$rm1,$"'d' 3v47_TLR5_Complex_Chains_A_D_Interactions_Only_SS.txt > 3v47_TLR5_Complex_Chains_A_D_Bonds_Only_SS.txt; fi
cut 3v47_TLR5_Complex_Chains_A_D_Bonds_Only_SS.txt -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > 3v47_TLR5_Complex_Chains_A_D_Bonds_A_Residues.txt

grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" 3v47_TLR5_Complex_Chains_B_C_Interactions.txt  | while read p; do echo $p; done > 3v47_TLR5_Complex_Chains_B_C_Interactions_Only_SS.txt
mapfile -t my_array < <( grep -n "^1. " 3v47_TLR5_Complex_Chains_B_C_Interactions_Only_SS.txt | cut -f1 -d: )
if [ ${#my_array[@]} == 3 ]; then rm1="${my_array[1]}"; rm2=$(echo "${my_array[2]}"-1 |bc); sed -e "$rm1,$rm2"'d' 3v47_TLR5_Complex_Chains_B_C_Interactions_Only_SS.txt > 3v47_TLR5_Complex_Chains_B_C_Bonds_Only_SS.txt; else rm1="${my_array[1]}"; sed -e "$rm1,$"'d' 3v47_TLR5_Complex_Chains_B_C_Interactions_Only_SS.txt > 3v47_TLR5_Complex_Chains_B_C_Bonds_Only_SS.txt; fi
cut 3v47_TLR5_Complex_Chains_B_C_Bonds_Only_SS.txt -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > 3v47_TLR5_Complex_Chains_B_C_Bonds_B_Residues.txt

grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" 3v47_TLR5_Complex_Chains_B_D_Interactions.txt  | while read p; do echo $p; done > 3v47_TLR5_Complex_Chains_B_D_Interactions_Only_SS.txt
mapfile -t my_array < <( grep -n "^1. " 3v47_TLR5_Complex_Chains_B_D_Interactions_Only_SS.txt | cut -f1 -d: )
if [ ${#my_array[@]} == 3 ]; then rm1="${my_array[1]}"; rm2=$(echo "${my_array[2]}"-1 |bc); sed -e "$rm1,$rm2"'d' 3v47_TLR5_Complex_Chains_B_D_Interactions_Only_SS.txt > 3v47_TLR5_Complex_Chains_B_D_Bonds_Only_SS.txt; else rm1="${my_array[1]}"; sed -e "$rm1,$"'d' 3v47_TLR5_Complex_Chains_B_D_Interactions_Only_SS.txt > 3v47_TLR5_Complex_Chains_B_D_Bonds_Only_SS.txt; fi
cut 3v47_TLR5_Complex_Chains_B_D_Bonds_Only_SS.txt -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > 3v47_TLR5_Complex_Chains_B_D_Bonds_B_Residues.txt

#Phantom script to reset Notepad++ $ for coloring
sed -e "$rm1,$w"'d'

#Combine the results into a per chain list (A or B) for bonds and interactions 
cat 3v47_TLR5_Complex_Chains_A_C_Bonds_A_Residues.txt 3v47_TLR5_Complex_Chains_A_D_Bonds_A_Residues.txt | sort > 3v47_TLR5_Complex_Bonds_Chain_A_Residues.txt
cat 3v47_TLR5_Complex_Chains_A_C_Interactions_A_Residues.txt 3v47_TLR5_Complex_Chains_A_D_Interactions_A_Residues.txt | sort > 3v47_TLR5_Complex_Interactions_Chain_A_Residues.txt


cat 3v47_TLR5_Complex_Chains_B_C_Bonds_B_Residues.txt 3v47_TLR5_Complex_Chains_B_D_Bonds_B_Residues.txt | sort > 3v47_TLR5_Complex_Bonds_Chain_B_Residues.txt
cat 3v47_TLR5_Complex_Chains_B_C_Interactions_B_Residues.txt 3v47_TLR5_Complex_Chains_B_D_Interactions_B_Residues.txt | sort > 3v47_TLR5_Complex_Interactions_Chain_B_Residues.txt

#2z7z-Ligand
#Create total interactions residues list for ligand to the reference
grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" 2z7x_TLR1_2_Complex_Ligand_Interactions.txt | while read p; do echo $p; done | cut -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq > 2z7x_TLR1_2_Complex_Ligand_Interactions_All_Residues.txt

#Separate the results into A and B chains
grep A 2z7x_TLR1_2_Complex_Ligand_Interactions_All_Residues.txt > 2z7x_TLR1_2_Complex_Ligand_Interactions_A_Residues.txt
grep B 2z7x_TLR1_2_Complex_Ligand_Interactions_All_Residues.txt > 2z7x_TLR1_2_Complex_Ligand_Interactions_B_Residues.txt


#Create total interactions residues list for each chain in the reference
grep -e "^[0-9]" -e "^ [0-9]" -e "^  [0-9]" 2z7x_TLR1_2_Complex_Ligand_Interactions.txt  | while read p; do echo $p; done > 2z7x_TLR1_2_Complex_Ligand_Interactions_All_Only_SS.txt
mapfile -t my_array < <( grep -n "^1. " 2z7x_TLR1_2_Complex_Ligand_Interactions_All_Only_SS.txt | cut -f1 -d: )
if [ ${#my_array[@]} == 3 ]; then rm1="${my_array[1]}"; rm2=$(echo "${my_array[2]}"-1 |bc); sed -e "$rm1,$rm2"'d' 2z7x_TLR1_2_Complex_Ligand_Interactions_All_Only_SS.txt > 2z7x_TLR1_2_Complex_Ligand_Bonds_Only_SS.txt; else rm1="${my_array[1]}"; sed -e "$rm1,$"'d' 2z7x_TLR1_2_Complex_Ligand_Interactions_All_Only_SS.txt > 2z7x_TLR1_2_Complex_Ligand_Bonds_Only_SS.txt; fi
cut 2z7x_TLR1_2_Complex_Ligand_Bonds_Only_SS.txt -f 5,6 -d ' ' | sed 's/ //g' | sort | uniq >  2z7x_TLR1_2_Complex_Ligand_Bonds_All_Residues.txt

#Phantom script to reset Notepad++ $ for coloring
sed -e "$rm1,$w"'d'

#Separate the results into A and B chains
grep A 2z7x_TLR1_2_Complex_Ligand_Bonds_All_Residues.txt > 2z7x_TLR1_2_Complex_Ligand_Bonds_A_Residues.txt
grep B 2z7x_TLR1_2_Complex_Ligand_Bonds_All_Residues.txt > 2z7x_TLR1_2_Complex_Ligand_Bonds_B_Residues.txt


#Download the top 30 pdb files from ClusPro and unzip them

#Use rename to add a unique identifier to the beginning of the file ffr
rename 's/model/TLR5_Homodimer_Any_Loc_Refine_Balanced_ClusPro_Model/' *.pdb

#Upload each of the top 10 structures to the Generate tool of PDBSum and input Email
#Once you get the email with the results, open it and move to the Prot-prot tab
#Once there click on the Interface Summary or Summaries you are interested on the left side
#Then scroll down to the "Residue interactions across interface" and open the "List of interactions" in a new tab
#Copy and paste all of this to a new txt file and save it with a unique Name
#Move everything to a new folder

#Make a uniquely named file to hold docking percentages and add column names
echo "TLR_Construct_Model" "Chain_A_Interactions_Percentage" "Chain_B_Interactions_Percentage" > 3v47_TLR5_Complex_Any_Loc_Docking_Percentages.txt

#Edit the draft_summarizing_pdbsum_interactions to your needs and copy over to terminal to run

#Copy the Total interactions to folder focused on A Chain bonds and rename them
rename 's/Interactions/Chain_A_Bonds/' *.txt

#Remove any B Chain and non-bonded data

#Make a uniquely named file to hold docking percentages and add column names
echo "TLR_Construct_Model" "Chain_A_Bonds_Percentage" > 3v47_TLR5_Complex_Any_Loc_Chain_A_Bond_Docking_Percentages.txt

#Edit the draft_chain_x_pdsum_bonds.sh to your needs and copy over to terminal to run 

#Copy the Total interactions to folder focused on B Chain bonds and rename them
rename 's/Interactions/Chain_B_Bonds/' *.txt

#Remove any A Chain and non-bonded data

#Make a uniquely named file to hold docking percentages and add column names
echo "TLR_Construct_Model" "Chain_B_Bonds_Percentage" > 3v47_TLR5_Complex_Any_Loc_Chain_B_Bond_Docking_Percentages.txt

#Edit the draft_chain_x_pdsum_bonds.sh to your needs and copy over to terminal to run 

#Combine all Docking Percentages Data into one Excel file and use to choose the best docking model (i.e. highest percentage of interactions and bonds)

##Molecular Dynamics

#Install GROMACS
#wget ftp://ftp.gromacs.org/gromacs/gromacs-2023.3.tar.gz
#tar xvf gromacs-2022.3.tar.gz
#cd gromacs-2022.3
#mkdir build
#cd build
#cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_INSTALL_PREFIX=/project/fsepru113/dbradshaw/gromacs-2023.3
#make
#make check
#make install

#Download default mdp files from compchems for reference
wget https://www.compchems.com/gromacs_protein_water/ions.mdp
wget https://www.compchems.com/gromacs_protein_water/minim.mdp
wget https://www.compchems.com/gromacs_protein_water/nvt.mdp
wget https://www.compchems.com/gromacs_protein_water/npt.mdp
wget https://www.compchems.com/gromacs_protein_water/md.mdp

#Start an interactive SLURM run
srun -N 1 --mem=144gb -c 72 -p short --pty $SHELL

#Activate gromacs
source /project/fsepru113/dbradshaw/gromacs-2023.3/bin/GMXRC

#Upload the molecular docking model you are testing
scp *.pdb david.bradshaw@ceres.scinet.usda.gov:/project/fsepru113/dbradshaw/reverse_vaccinology/Molecular_Dynamics/Loc_Indpdt_TLR2

#Open draft_md_pre_production_run.sh and adjust as needed then copy over to HPC and run it
nano Loc_Indpdt_TLR5_Complex_ClusPro_Model_4_md_pre_production_run.sh

#Delete water and ligands from pdb file
grep -v HETATM TLR2_Any_Loc_Refined_Balanced_ClusPro_model.000.00.pdb > TLR2_Any_clean.pdb

#Generate a topology file - < 1 min
gmx pdb2gmx -f TLR2_Any_clean.pdb -o TLR2_Any_clean.gro -water spce -ignh
#Choose option 15 for OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)

#Define a box - < 1 min
gmx editconf -f TLR2_Any_clean.gro -o TLR2_Any_box.gro -c -d 1.0 -bt cubic

#Solvate the system - < 1 min
gmx solvate -cp TLR2_Any_box.gro -cs spc216.gro -o TLR2_Any_solv.gro -p topol.top

##Neutralize the system

#Assemble the ions.tpr file  - < 1 min
gmx grompp -f ions.mdp -c TLR2_Any_solv.gro -p topol.top -o ions.tpr

#Genion module  - < 1 min
gmx genion -s ions.tpr -o TLR2_Any_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
#Choose option 13 - SOL


##Energy minimization

#Generate tpr file - < 1 min
gmx grompp -f minim.mdp -c TLR2_Any_ions.gro -p topol.top -o em.tpr

#Run simulation - < 3 hrs
gmx mdrun -v -deffnm em

#Generate a plot of potential over time
gmx energy -f em.edr -o potential.xvg


##NVT equilibration

#Generate a tpr file - < 1 min
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

#Run simulation 
gmx mdrun -v -deffnm nvt

#Generate a graphic for GROMACS energy - < 1 min
printf "16 0" | gmx energy -f nvt.edr -o nvt_temperature.xvg
#Choose 16 for Temperature

#Generate a graphic for Root mean square deviation (RMSD) - < 1 min
printf "4 4" | gmx rms -f nvt.trr -s nvt.tpr -o nvt_rmsd.xvg
#Choose 4 to select backbone group for least squares fit and 4 to choose backbone group for RMSD calculation


##NPT equilibration


#Generate a tpr file - < 1 min
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

#Run simulation 
gmx mdrun -v -deffnm npt

#Generate a graphic for GROMACS energy - < 1 min
printf "18 0" | gmx energy -f npt.edr -o npt_pressure.xvg
#Choose 18 for Pressure

#Generate a graphic for RMSD - < 1 min
printf "4 4" | gmx rms -f npt.trr -s npt.tpr -o npt_rmsd.xvg
#Choose 4 to select backbone group for least squares fit and 4 to choose backbone group for RMSD calculation


##Production Run

#Use a interactive srun script to determine how long analysis will take, change as you troubleshoot the resources
#srun -N 1 --mem=72gb -c 36 -p short --pty $SHELL
srun -N 1 --mem=144gb -c 72 -p short --pty $SHELL

#Activate gromacs
source /project/fsepru113/dbradshaw/gromacs-2023.3/bin/GMXRC

#Use nano to change the nsteps in md.mdp as needed (starting with 10 ns)
nano md.mdp

#Generate a tpr file - < 1 min
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

#Run a test simulation
gmx mdrun -v -deffnm md

#Pay attention to the screen, the predicted finishing time will quickly flash across the screen
#Example - running in srun Loc_Indpdt_TLR2 on 12/6/23 at 13:25 would complete on 12/23/23 at 2:30 using 72 cores and 144gb i.e. ~18 days (long partition is 21 days, so thats the target)

#Use Ctrl+C to kill the process
#Remove the resulting uncomplete files (use ls -lah to double check time if unsure)
rm md.edr
rm md.log
rm md.xtc

#Either adjust your resources and repeat the above or if satisfied adjust draft_md_production_run.sh to match your needs and copy over to new nano .sh file
nano Loc_Indpdt_TLR2_md_prod_run.sh
sbatch Loc_Indpdt_TLR2_md_prod_run.sh --no-requeue

##Restarting the Run

#-cpi
#-v is to make it verbose
#-cpi is to restart from the latest checkpoint
#-cpt is to adjust the frequency of the checkpoint output in minutes
#Adjust the draft_md_production_run.sh to have the following
gmx mdrun -v -deffnm md -cpi md.cpt

sbatch Loc_Indpdt_TLR2_md_prod_run.sh --no-requeue


#From https://www.compchems.com/extend-or-continue-a-gromacs-simulation/#continue-a-simulation
#Note that if you restart different simulations from the same checkpoint you will find that the continuations will diverge
#from each other. This is due to the limited precision of computers at our disposal. In principle, you should be able to
#exactly reproduce the same results with an optimal computer having unlimited precision.
#Despite this, different trajectories are all equally valid and none of them is better than the others. So don’t worry about this.

##Generate graphics

#Generate a graphic for RMSD - < 5 min (10 ns sim)
#https://www.compchems.com/what-is-the-rmsd-and-how-to-compute-it-with-gromacs/
#-s is the reference geometry
#-f is the trajectory file
#-tu ns will change the unites to ns instead of ps
printf "4 4" | gmx rms -f md.xtc -s md.tpr -o md_rmsd.xvg
#Choose 4 to select backbone group for least squares fit and 4 to choose backbone group for RMSD calculation

#Generate a graphic for radius of gyration - < 3 min (10 ns sim)
#https://tutorials.gromacs.org/docs/md-intro-tutorial.html
printf "1" | gmx gyrate -f md.xtc -s md.tpr -o md_gyrate.xvg
#Choose 1  to select protein group

#Generate a graphic for RMSF (Root Mean Square Fluctuation) - < 3 min (10 ns sim)
#https://www.compchems.com/how-to-compute-the-rmsf-using-gromacs/#what-is-the-rmsf
#From https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html - An area of the structure with high RMSF values frequently diverges from the average, indicating high mobility. When RMSF analysis is carried out on proteins, it is typically restricted to backbone or alpha-carbon atoms; these are more characteristic of conformational changes than the more flexible side-chains.
#From https://gromacs.bioexcel.eu/t/rmsf-analysis-of-protein-ligand-md-simulation/4919/5 - I would also recommend performing RMSF analysis only on the Cα or backbone groups, not the whole protein. You will get a lot of irrelevant motions (like sidechains) if you consider everything.
# -f is the trajectory file
# -s is the reference structure
# -res is to compute the RMSF for each Residue
# -b and -e allow you to choose a specific time frame in ps
printf "3" | gmx rmsf -f md.xtc -s md.tpr -o md_rmsf.xvg -res
#Choose 3 to select the C-alpha group

###Codon Optimization
#JCAT 
#https://jcat.de
#Made RE pair choices based upon the multiple cloning site (MCS) of pET-30a vector 
#Copy the protein sequence of the construct into the input area
#Choose Protein option
#Select Escherichia coli (strain K12) in teh dropdown menu below
#Click "Avoid Clevage Sites of Restriction Enzymes" and choose the enzymes you are testing
#Tested the following pairs (3' and 5'): XhoI & BamHI, EcoRI & BamHI, and HindIII & BamHI
#Submit and choose the RE pair with highest CAI (preferably 1.0) or test different pairs based upon availability in the the MCS and support from the literature
#Save the output as a pdf
#Record the GC content for the E. coli K12 Strain and the corrected sequence and RE pair you have chosen
#Copy and paste the sequence to Excel, use Text to Column to remove the sequence from line numbers, remove the white space with Find and Replace, and use then TextJoin formula to concatenate the sequences


###In silico cloning

##Input Files
#Download the .dna from SnapGene via File -> Open Files
#https://www.snapgene.com/plasmids/pet_and_duet_vectors_(novagen)/pET-30a(%2B)
#Copy and paste the codon optimized sequence into Notepad++ and add the restriction enzyme sites (XhoI to N terminal and BamHI to C terminal) to the ends and save the file

##Genious
#Go to Add -> Import Files -> Choose construct fasta file and vector .dna file
#Check both options for imporing SnapGene Sequence file (Also import enzme sets associated with <vector> & Annotate teh restriction enzymes that were displayed in SnapGene)
#Click both the construct and vector and go to Cloning -> Restriction Cloning
#Backbone is pET-30a(+) with the BamHI RE site as the 5' end and the XhoI RE site as the 3' ends
#Insert is the construct with the XhoI RE site as the 5' end and the BamHI RE site as teh 5' end
#Click Generate Constructs
#Click the resulting file and then go to Export -> To Multiple Files, Choose folder, save as a .geneious file

##SnapGene Viewer
#Go to Open -> find .geneious file and load it in
#Change the color and name of Construct Insert by double clicking it 
#While selecting the Construct Insert, right click the DNA Sequence and select Set DNA Color -> Do both bands
#Right click enzymes used and then select Highlight Enzyme Site
#Double click name in center to change it to something like  "pET-30a(+) - Codon Optimized Construct Sequence - XhoI BamHI"
#Use drop down next to Save to Save as the entire file as SnapGene DNA file
#Use drop down next to Save to Export Map as a TIFF with 300 dpi


###In Silico Immune Trials
#C-ImmSim
#https://kraken.iac.rm.cnr.it/C-IMMSIM/
#https://wwwold.iac.rm.cnr.it/~filippo/c-immsim/index.html
#https://wwwold.iac.rm.cnr.it/~filippo/c-immsim/the-parameters.html
#https://wwwold.iac.rm.cnr.it/~filippo/c-immsim/ewExternalFiles/cimmsim-description-examples.pdf
#Keep Random Seed at 12345 and Simulation Voume at 10
#Each step is 8 hours thus for a 12 week trial need to change Simulations Steps to 252 steps (12 weeks x 7 = 84 days; 84 days * 24 = 2016 hours; 2016 hours / 8 = 252 steps)
#Change MHC molecules to match what was used in NetMHCpan and NetMCHIIpan in pairs of doubles
#for injections on Day 0 and 21 do the following:
#Injection N. 1 = Time Steop of Injection = 1 (8 hours, Day 0); What to inject: to vaccine (no LPS) & Adjuvant to 100; Num Ag to inject to 100; Copy and paste in protein sequence of Construct
#Injection N. 2 = Time Steop of Injection = 63 (21 x 24 / 8); What to inject: to vaccine (no LPS) & Adjuvant to 100; Num Ag to inject to 100; Copy and paste in protein sequence of Construct
#Click Submit Job
#Once finished click on Output, then in the output click Download PDF Report, and save as unique file


####ChimeraX####
#https://www.cgl.ucsf.edu/chimerax/
#https://www.youtube.com/playlist?list=PL4eF1KHNgDfIYSKCS3_S0PTRYtYTV9Myi
#https://www.cgl.ucsf.edu/chimerax/docs/credits.html

#Change information about chains
#https://www.cgl.ucsf.edu/chimerax/docs/user/commands/changechains.html

###Any Loc Construct - Original

##Add colors to construct based upon region
#Adjuvant
color #1/?:1-495 blue

#CTL Region
color #1/?:496-605 green

#Th Region
color #1/?:606-785 orange

#LBL Region
color #1/?:786-936 red

#Save image
#https://www.cgl.ucsf.edu/chimerax/docs/user/commands/save.html
#https://mail.cgl.ucsf.edu/mailman/archives/list/chimerax-users@cgl.ucsf.edu/message/LDM5SPK5HMV67QMEI6ZCKVCR7T7VAKXO/
save "I:\Bearson\David_Bradshaw\Manuscripts\UK1_Reverse_Vaccinology\Figures\Figure_3a_Original_Phyre_3D_Structure_V2.tif" width 2100 height 1200 transparentBackground true
#Target Size = Height x Width = 2" x 3.5"
#2 x 2 x 300 dpi = 1200 pixels height
#3.5 x 2 x 300 dpi = 2100 pixels width

#Save ChimeraX Session
save "I:\Bearson\David_Bradshaw\Manuscripts\UK1_Reverse_Vaccinology\Figures\Figure_3a_Original_Phyre_3D_Structure_V2.cxs"

###Any Loc Construct - Refined

##Add colors to construct based upon region
#Adjuvant
color #1/?:1-495 blue

#CTL Region
color #1/?:496-605 green

#Th Region
color #1/?:606-785 orange

#LBL Region
color #1/?:786-936 red

#Save image
#https://www.cgl.ucsf.edu/chimerax/docs/user/commands/save.html
save "I:\Bearson\David_Bradshaw\Manuscripts\UK1_Reverse_Vaccinology\Figures\Figure_3b_Refined_Model_4_3D_Structure_V2.tif" width 2100 height 1200 transparentBackground true
#Target Size = Height x Width = 2" x 3.5"
#2 x 2 x 300 dpi = 1200 pixels height
#3.5 x 2 x 300 dpi = 2100 pixels width

#Save ChimeraX Session
save "I:\Bearson\David_Bradshaw\Manuscripts\UK1_Reverse_Vaccinology\Figures\Figure_3b_Refined_Model_4_3D_Structure_V2.cxs"

###Any Loc Construct - TLR 1/2


##Add colors to construct based upon region
#Adjuvant
color #1/?:1-495 blue

#CTL Region
color #1/?:496-605 green

#Th Region
color #1/?:606-785 orange

#LBL Region
color #1/?:786-936 red

#Add colors to TLRs by Chain name

#TLR-1
color /a cyan

#TLR-2
color /b yellow

#Save image
#https://www.cgl.ucsf.edu/chimerax/docs/user/commands/save.html
save "I:\Bearson\David_Bradshaw\Manuscripts\UK1_Reverse_Vaccinology\Figures\Figure_4a_TLR1_2_Heterodimer_Any_Loc_Refined_Balanced_ClusPro_Model_2_V2.tif" width 2100 height 1200 transparentBackground true
#Target Size = Height x Width = 2" x 3.5"
#2 x 2 x 300 dpi = 1200 pixels height
#3.5 x 2 x 300 dpi = 2100 pixels width

#Save ChimeraX Session
save "I:\Bearson\David_Bradshaw\Manuscripts\UK1_Reverse_Vaccinology\Figures\Figure_4a_TLR1_2_Heterodimer_Any_Loc_Refined_Balanced_ClusPro_Model_2_V2.cxs"

###Any Loc Construct - TLR 5

##Add colors to construct based upon region
#Adjuvant
color #1/?:1-495 blue

#CTL Region
color #1/?:496-605 green

#Th Region
color #1/?:606-785 orange

#LBL Region
color #1/?:786-936 red

#Add colors to TLRs by Chain name

#TLR-5A
color /a pink

#TLR-5B
color /b violet

#Save image
#https://www.cgl.ucsf.edu/chimerax/docs/user/commands/save.html
save "I:\Bearson\David_Bradshaw\Manuscripts\UK1_Reverse_Vaccinology\Figures\Figure_4b_TLR5_Homodimer_Any_Loc_Refined_Balanced_ClusPro_Model_6_V2.tif" width 2100 height 1200 transparentBackground true
#Target Size = Height x Width = 2" x 3.5"
#2 x 2 x 300 dpi = 1200 pixels height
#3.5 x 2 x 300 dpi = 2100 pixels width

#Save ChimeraX Session
save "I:\Bearson\David_Bradshaw\Manuscripts\UK1_Reverse_Vaccinology\Figures\Figure_4b_TLR5_Homodimer_Any_Loc_Refined_Balanced_ClusPro_Model_6_V2.cxs"

##Adjust DPI of an image
#https://guides.lib.umich.edu/c.php?g=282942&p=1888164
#Right-click -> Open with -> GIMP or similar software
#Crop as needed -> Default mouse tool upon opening -> Highlight area want to keep -> press Enter
#GIMP -> Image -> Print Size -> Adjust DPI to desired level, size will decrease if increasing dpi
#GIMP -> File -> Export As -> ""_300dpi.png -> Export with default settings

##Check dpi of a PNG
#https://answers.microsoft.com/en-us/windows/forum/all/how-to-check-dpi-on-png/3e18a21c-16ae-464e-b80d-a7d73c2e121c
#Right-click -> Open With -> Paint 
#Paint -> File -> Image Properties