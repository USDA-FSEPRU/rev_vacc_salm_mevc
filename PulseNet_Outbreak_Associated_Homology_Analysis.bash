###Outbreak Comparisons

#Make directories for each of the sets of data
mkdir Set_1_Enteritidis
mkdir Set_2_Newport_Typhimurium
mkdir Set_3_Other_GT_1000
mkdir Set_4_LT_1000_GT_250
mkdir Set_5_LT_250

#Send the list of SRRs to the HPC and mv them to appropriate folders
mv PulseNet_Newport_Typhimurium_SRRs.txt Set_2_Newport_Typhimurium/
mv PulseNet_other_GT_1000_SRRs.txt Set_3_Other_GT_1000/
mv PulseNet_LT_1000_GT_250_SRRs.txt Set_4_LT_1000_GT_250/
mv PulseNet_LT_250_SRRs.txt Set_5_LT_250/


###Set 1 Enteritidis - 6245

##Step 1 SRRs to SRAs

#Create a srun
srun -N 1 --mem=48gb -c 24 -p short --pty $SHELL

#Get all the reads using step_1_draft_srr_to_sras.sh n= 1509
nano  Step_1_PulseNet_Enteritidis_SRRs.sh
sbatch  Step_1_PulseNet_Enteritidis_SRRs.sh

#Remove empty directories
rmdir SRR*

#Make a list of the sra files that got dounloaded
ls sras/ | cut -d "." -f 1 > Rd1_sras.txt
wc -l Rd1_sras.txt
#6211

#Make a missing summarized sras file and record the number
grep -Fxv -f Rd1_sras.txt PulseNet_Enteritidis_SRRs.txt> missing_sras_post_downloading_Rd1.txt
wc -l missing_sras_post_downloading_Rd1.txt
#34

#Load the sra toolkit
module load sratoolkit

#Grab the missing files "locally" to see if any erros occur
prefetch --option-file missing_sras_post_downloading_Rd1.txt --verbose

#Move sras to the sra folder
cat missing_sras_post_downloading_Rd1.txt | while read srr; do mv /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/"$srr"/* sras/; done
rmdir SRR*
ls sras/ | wc -l
#6244

#SRR164787
#No data for this SRR
#6244 is the new expected total for Enteritidis

#Make a new upated file based upon sras directory
ls sras/ | cut -d "." -f 1 > PulseNet_Enteritidis_SRRs_updated.txt
wc -l PulseNet_Enteritidis_SRRs_updated.txt
#6244

#Make a new directory for the reads and folders in reads folder to separate results
mkdir reads
mkdir reads/raw
mkdir reads/fastp
mkdir reads/unpaired

#Change all sras to paired reads using step_2_draft_sra_to_reads.sh
nano Step_2_PulseNet_Enteritidis_sras_to_reads.sh
sbatch Step_2_PulseNet_Enteritidis_sras_to_reads.sh

#Make a list of unique ids of the sras that got split into reads (and double check you got right number)
ls reads/raw/*_1.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R1_reads.txt
ls reads/raw/*_2.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R2_reads.txt

wc -l Rd1_R1_reads.txt
#6240

wc -l Rd1_R2_reads.txt
#6244

#Make a missing summarized sras file and record the number
grep -Fxv -f Rd1_R1_reads.txt PulseNet_Enteritidis_SRRs_updated.txt> missing_reads_post_Rd1.txt
wc -l missing_reads_post_Rd1.txt
#4


#Check these SRRs in the SRAs
#All four are single read runs, thus make a list of them
grep -Fxv -f Rd1_R1_reads.txt Rd1_R2_reads.txt> single_end_SRRs.txt
wc single_end_SRRs.txt
#4

#Remove single end reads
cat single_end_SRRs.txt | while read sra; do rm reads/raw/"$sra"*; done

#Make a new upated file based upon reads directory with paired end reads
ls reads/raw/*_1.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > PulseNet_Enteritidis_SRRs_paired_end.txt
wc -l PulseNet_Enteritidis_SRRs_paired_end.txt
#6240

#Run fastp on all reads without assemblies using step_3_draft_fastp_raw_reads.sh
nano Step_3_PulseNet_Enteritidis_fastp_raw_reads.sh
sbatch Step_3_PulseNet_Enteritidis_fastp_raw_reads.sh

#Make a list of unique ids of the reads that fastp worked on
ls reads/fastp/*  | cut -d "/" -f 3 | cut -d "." -f 1 | uniq > Rd1_fastp.txt
wc -l Rd1_fastp.txt
#6232/6240

#Make a missing fastp reads file and record the number
grep -Fxv -f Rd1_fastp.txt PulseNet_Enteritidis_SRRs_paired_end.txt> missing_fastp_post_Rd1.txt
wc -l missing_fastp_post_Rd1.txt
#8

#Activate an environment with fastp
source /project/fsepru113/dbradshaw/miniconda3/bin/activate polish

#Run Fastp 
cat missing_fastp_post_Rd1.txt | while read read_id; do fastp --in1 reads/raw/"$read_id"_1.fastq --in2 reads/raw/"$read_id"_2.fastq --out1 reads/fastp/"$read_id".fastp_1.fastq --out2 reads/fastp/"$read_id".fastp_2.fastq --unpaired1 reads/unpaired/"$read_id".fastp.unpaired_1.fastq --unpaired2 reads/unpaired/"$read_id"_fastp.unpaired_2.fastq; done

#Check that you have right number
ls reads/fastp/*_1.fastq | wc -l 
#6240/6240

ls reads/unpaired/*_1.fastq | wc -l 
#6240/6240

###SKESA

#Make a directory
mkdir skesa_assembly/

#Adjust step_4a_draft_skesa_assembly_reads_folder.sh to your situation
nano Step_4_PulseNet_Enteritidis_skesa_assembly.sh
sbatch Step_4_PulseNet_Enteritidis_skesa_assembly.sh


#Check for correctly made assemblies, do by size becasue will output even if it fails
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#6226/6240


#Check for skesa assemblies that were not made correctly
find skesa_assembly/*.fna -type f -size -1c -exec mv "{}" redo_assemblies/ \;

#Determine how many contigs were moved to that folder
ls redo_assemblies/ | wc -l
#14/6240

#Make a list of unique ids of the reads that fastp worked on
ls redo_assemblies/ | cut -d "." -f 1 > missing_assemblies_post_Rd1.txt
wc -l missing_assemblies_post_Rd1.txt
#14

#Move assemblies that did work to final_assemblies
mv skesa_assembly/*.fna final_assemblies
ls final_assemblies | wc -l
#6226/6240

#Activate an environment with skesa
source /project/fsepru113/dbradshaw/miniconda3/bin/activate skesa

#Run Fastp 
cat missing_assemblies_post_Rd1.txt | while read read_id; do skesa --reads reads/fastp/"$read_id".fastp_1.fastq,reads/fastp/"$read_id".fastp_2.fastq --cores 4 --memory 48 > skesa_assembly/"$read_id".skesa.fna; done

#Check that you have right number
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#6240/6240

#Create a final txt to hold all SRRs that made it through the process
cp PulseNet_Enteritidis_SRRs_paired_end.txt PulseNet_Enteritidis_SRRs_final.txt

#Move the files to its appropriate Set
mv (*.txt,*.sh) Set_1_Enteritidis/

#Remove the sras and reads but keep foldes
rm sras/*
rm reads/raw/*
rm reads/fastp/*
rm reads/unpaired/*

#Move the assemblies from temp folder to final folder
mv skesa_assembly/* final_assemblies

###Set 2 Newport_Typhimurium - 7767

##Step 1 SRRs to SRAs

#Move over the txt file to main directory
mv Set_2_Newport_Typhimurium/PulseNet_Newport_Typhimurium_SRRs.txt .

#Create a srun
srun -N 1 --mem=48gb -c 24 -p short --pty $SHELL

#Get all the reads using step_1_draft_srr_to_sras.sh n= 1509
nano  Step_1_PulseNet_Newport_Typhimurium_SRRs.sh
sbatch  Step_1_PulseNet_Newport_Typhimurium_SRRs.sh

#Remove empty directories
rmdir SRR*

#Make a list of the sra files that got dounloaded
ls sras/ | cut -d "." -f 1 > Rd1_sras.txt
wc -l Rd1_sras.txt
#7726/7767

#Make a missing summarized sras file and record the number
grep -Fxv -f Rd1_sras.txt PulseNet_Newport_Typhimurium_SRRs.txt> missing_sras_post_downloading_Rd1.txt
wc -l missing_sras_post_downloading_Rd1.txt
#41

#Load the sra toolkit
module load sratoolkit

#Grab the missing files "locally" to see if any erros occur
prefetch --option-file missing_sras_post_downloading_Rd1.txt --verbose

#Move sras to the sra folder
cat missing_sras_post_downloading_Rd1.txt | while read srr; do mv /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/"$srr"/* sras/; done
rmdir SRR*
ls sras/ | wc -l
#7766/7767

#SRR164485
#No data for this SRR
#7766 is the new expected total for Newport_Typhimurium

#Make a new upated file based upon sras directory
ls sras/ | cut -d "." -f 1 > PulseNet_Newport_Typhimurium_SRRs_updated.txt
wc -l PulseNet_Newport_Typhimurium_SRRs_updated.txt
#7766

#Change all sras to paired reads using step_2_draft_sra_to_reads.sh
nano Step_2_PulseNet_Newport_Typhimurium_sras_to_reads.sh
sbatch Step_2_PulseNet_Newport_Typhimurium_sras_to_reads.sh

#Make a list of unique ids of the sras that got split into reads (and double check you got right number)
ls reads/raw/*_1.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R1_reads.txt
ls reads/raw/*_2.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R2_reads.txt

wc -l Rd1_R1_reads.txt
#7766

wc -l Rd1_R2_reads.txt
#7766

#No missing reads

#Run fastp on all reads without assemblies using step_3_draft_fastp_raw_reads.sh
nano Step_3_PulseNet_Newport_Typhimurium_fastp_raw_reads.sh
sbatch Step_3_PulseNet_Newport_Typhimurium_fastp_raw_reads.sh

#Make a list of unique ids of the reads that fastp worked on
ls reads/fastp/*  | cut -d "/" -f 3 | cut -d "." -f 1 | uniq > Rd1_fastp.txt
wc -l Rd1_fastp.txt
#7763/7766

#Make a missing fastp reads file and record the number
grep -Fxv -f Rd1_fastp.txt PulseNet_Newport_Typhimurium_SRRs_updated.txt> missing_fastp_post_Rd1.txt
wc -l missing_fastp_post_Rd1.txt
#3

#Activate an environment with fastp
source /project/fsepru113/dbradshaw/miniconda3/bin/activate polish

#Run Fastp 
cat missing_fastp_post_Rd1.txt | while read read_id; do fastp --in1 reads/raw/"$read_id"_1.fastq --in2 reads/raw/"$read_id"_2.fastq --out1 reads/fastp/"$read_id".fastp_1.fastq --out2 reads/fastp/"$read_id".fastp_2.fastq --unpaired1 reads/unpaired/"$read_id".fastp.unpaired_1.fastq --unpaired2 reads/unpaired/"$read_id"_fastp.unpaired_2.fastq; done

#Check that you have right number
ls reads/fastp/*_1.fastq | wc -l 
#7766/7766

ls reads/unpaired/*_1.fastq | wc -l 
#7766/7766

ls reads/fastp/*_2.fastq | wc -l 
#7766/7766

ls reads/unpaired/*_2.fastq | wc -l 
#7766/7766


###SKESA

#Make a directory
mkdir skesa_assembly/

#Adjust step_4a_draft_skesa_assembly_reads_folder.sh to your situation
nano Step_4_PulseNet_Newport_Typhimurium_skesa_assembly.sh
sbatch Step_4_PulseNet_Newport_Typhimurium_skesa_assembly.sh


#Check for correctly made assemblies, do by size becasue will output even if it fails
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#7762/7766

#Remove assemblies that are less than 1 byte
find skesa_assembly/*.fna -type f -size -1c -exec rm "{}" \;

#Get a list of assemblies that worked the first time
ls skesa_assembly/ | cut -d "." -f 1 > assemblies_Rd1.txt
wc -l assemblies_Rd1.txt
#7762

#Make a list of unique ids of the reads that fastp worked on
grep -Fxv -f assemblies_Rd1.txt PulseNet_Newport_Typhimurium_SRRs_updated.txt> missing_assemblies_post_Rd1.txt
wc -l missing_assemblies_post_Rd1.txt
#4

#Activate an environment with skesa
source /project/fsepru113/dbradshaw/miniconda3/bin/activate skesa

#Run Fastp 
cat missing_assemblies_post_Rd1.txt | while read read_id; do skesa --reads reads/fastp/"$read_id".fastp_1.fastq,reads/fastp/"$read_id".fastp_2.fastq --cores 4 --memory 48 > skesa_assembly/"$read_id".skesa.fna; done

#Check that you have right number
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#7766/7766

#Create a final txt to hold all SRRs that made it through the process
cp PulseNet_Newport_Typhimurium_SRRs_updated.txt PulseNet_Newport_Typhimurium_SRRs_final.txt

#Move the files to its appropriate Set
mv (*.txt,*.sh) Set_2_Newport_Typhimurium/

#Move the assemblies from temp folder to final folder
mv skesa_assembly/* final_assemblies
find final_assemblies/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#14006

#Remove the sras and reads but keep foldes
rm sras/*
rm reads/raw/*
rm reads/fastp/*
rm reads/unpaired/*


###Set 3 other_GT_1000 - 6927
mv 
##Step 1 SRRs to SRAs

#Move over the txt file to main directory
mv Set_2_other_GT_1000/PulseNet_other_GT_1000_SRRs.txt .

#Create a srun
srun -N 1 --mem=48gb -c 24 -p short --pty $SHELL

#Get all the reads using step_1_draft_srr_to_sras.sh n= 1509
nano  Step_1_PulseNet_other_GT_1000_SRRs.sh
sbatch  Step_1_PulseNet_other_GT_1000_SRRs.sh

#Remove empty directories
rmdir SRR*

#Make a list of the sra files that got dounloaded
ls sras/ | cut -d "." -f 1 > Rd1_sras.txt
wc -l Rd1_sras.txt
#6899/6927

#Make a missing summarized sras file and record the number
grep -Fxv -f Rd1_sras.txt PulseNet_other_GT_1000_SRRs.txt> missing_sras_post_downloading_Rd1.txt
wc -l missing_sras_post_downloading_Rd1.txt
#28

#Load the sra toolkit
module load sratoolkit

#Grab the missing files "locally" to see if any erros occur
prefetch --option-file missing_sras_post_downloading_Rd1.txt --verbose

#Move sras to the sra folder
cat missing_sras_post_downloading_Rd1.txt | while read srr; do mv /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/"$srr"/* sras/; done
rmdir SRR*
ls sras/ | wc -l
#6927/6927

#Change all sras to paired reads using step_2_draft_sra_to_reads.sh
nano Step_2_PulseNet_other_GT_1000_sras_to_reads.sh
sbatch Step_2_PulseNet_other_GT_1000_sras_to_reads.sh

#Make a list of unique ids of the sras that got split into reads (and double check you got right number)
ls reads/raw/*_1.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R1_reads.txt
ls reads/raw/*_2.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R2_reads.txt

wc -l Rd1_R1_reads.txt
#6927

wc -l Rd1_R2_reads.txt
#6927

#No missing reads

#Run fastp on all reads without assemblies using step_3_draft_fastp_raw_reads.sh
nano Step_3_PulseNet_other_GT_1000_fastp_raw_reads.sh
sbatch Step_3_PulseNet_other_GT_1000_fastp_raw_reads.sh

#Make a list of unique ids of the reads that fastp worked on
ls reads/fastp/*  | cut -d "/" -f 3 | cut -d "." -f 1 | uniq > Rd1_fastp.txt
wc -l Rd1_fastp.txt
#6914/6927

#Make a missing fastp reads file and record the number
grep -Fxv -f Rd1_fastp.txt PulseNet_other_GT_1000_SRRs.txt> missing_fastp_post_Rd1.txt
wc -l missing_fastp_post_Rd1.txt
#13

#Activate an environment with fastp
source /project/fsepru113/dbradshaw/miniconda3/bin/activate polish

#Run Fastp 
cat missing_fastp_post_Rd1.txt | while read read_id; do fastp --in1 reads/raw/"$read_id"_1.fastq --in2 reads/raw/"$read_id"_2.fastq --out1 reads/fastp/"$read_id".fastp_1.fastq --out2 reads/fastp/"$read_id".fastp_2.fastq --unpaired1 reads/unpaired/"$read_id".fastp.unpaired_1.fastq --unpaired2 reads/unpaired/"$read_id"_fastp.unpaired_2.fastq; done

#Check that you have right number
ls reads/fastp/*_1.fastq | wc -l 
#6927/6927

ls reads/unpaired/*_1.fastq | wc -l 
#6927/6927

ls reads/fastp/*_2.fastq | wc -l 
#6927/6927

ls reads/unpaired/*_2.fastq | wc -l 
#6927/6927


###SKESA

#Adjust step_4a_draft_skesa_assembly_reads_folder.sh to your situation
nano Step_4_PulseNet_other_GT_1000_skesa_assembly.sh
sbatch Step_4_PulseNet_other_GT_1000_skesa_assembly.sh


#Check for correctly made assemblies, do by size becasue will output even if it fails
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#6925/6927

#Remove assemblies that are less than 1 byte
find skesa_assembly/*.fna -type f -size -1c -exec rm "{}" \;

#Get a list of assemblies that worked the first time
ls skesa_assembly/ | cut -d "." -f 1 > assemblies_Rd1.txt
wc -l assemblies_Rd1.txt
#6925

#Make a list of unique ids of the reads that fastp worked on
grep -Fxv -f assemblies_Rd1.txt PulseNet_other_GT_1000_SRRs.txt> missing_assemblies_post_Rd1.txt
wc -l missing_assemblies_post_Rd1.txt
#4

#Activate an environment with skesa
source /project/fsepru113/dbradshaw/miniconda3/bin/activate skesa

#Run Fastp 
cat missing_assemblies_post_Rd1.txt | while read read_id; do skesa --reads reads/fastp/"$read_id".fastp_1.fastq,reads/fastp/"$read_id".fastp_2.fastq --cores 4 --memory 48 > skesa_assembly/"$read_id".skesa.fna; done

#Check that you have right number
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#6927/6927

#Create a final txt to hold all SRRs that made it through the process
cp PulseNet_other_GT_1000_SRRs.txt PulseNet_other_GT_1000_SRRs_final.txt

#Move the files to its appropriate Set
mv *.txt *.sh Set_3_Other_GT_1000/

#Move the assemblies from temp folder to final folder
mv skesa_assembly/* final_assemblies
find final_assemblies/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#20,933

#Remove the sras and reads but keep foldes
rm sras/*
rm reads/raw/*
rm reads/fastp/*
rm reads/unpaired/*

###Set 4 LT_1000_GT_250 - 5763

##Step 1 SRRs to SRAs

#Move over the txt file to main directory
mv Set_4_LT_1000_GT_250/PulseNet_LT_1000_GT_250_SRRs.txt .

#Create a srun
srun -N 1 --mem=48gb -c 24 -p short --pty $SHELL

#Get all the reads using step_1_draft_srr_to_sras.sh n= 1509
nano  Step_1_PulseNet_LT_1000_GT_250_SRRs.sh
sbatch  Step_1_PulseNet_LT_1000_GT_250_SRRs.sh

#Remove empty directories
rmdir SRR*

#Make a list of the sra files that got dounloaded
ls sras/ | cut -d "." -f 1 > Rd1_sras.txt
wc -l Rd1_sras.txt
#5553/5763

#Make a missing summarized sras file and record the number
grep -Fxv -f Rd1_sras.txt PulseNet_LT_1000_GT_250_SRRs.txt> missing_sras_post_downloading_Rd1.txt
wc -l missing_sras_post_downloading_Rd1.txt
#210

#Load the sra toolkit
module load sratoolkit

#Grab the missing files "locally" to see if any erros occur
prefetch --option-file missing_sras_post_downloading_Rd1.txt --verbose

#Move sras to the sra folder
cat missing_sras_post_downloading_Rd1.txt | while read srr; do mv /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/"$srr"/* sras/; done
rmdir SRR*
ls sras/ | wc -l
#5763/5763

#Change all sras to paired reads using step_2_draft_sra_to_reads.sh
nano Step_2_PulseNet_LT_1000_GT_250_sras_to_reads.sh
sbatch Step_2_PulseNet_LT_1000_GT_250_sras_to_reads.sh

#Make a list of unique ids of the sras that got split into reads (and double check you got right number)
ls reads/raw/*_1.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R1_reads.txt
ls reads/raw/*_2.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R2_reads.txt

wc -l Rd1_R1_reads.txt
#5763

wc -l Rd1_R2_reads.txt
#5763

#No missing reads

#Run fastp on all reads without assemblies using step_3_draft_fastp_raw_reads.sh
nano Step_3_PulseNet_LT_1000_GT_250_fastp_raw_reads.sh
sbatch Step_3_PulseNet_LT_1000_GT_250_fastp_raw_reads.sh

#Make a list of unique ids of the reads that fastp worked on
ls reads/fastp/*  | cut -d "/" -f 3 | cut -d "." -f 1 | uniq > Rd1_fastp.txt
wc -l Rd1_fastp.txt
#576/5763

#Make a missing fastp reads file and record the number
grep -Fxv -f Rd1_fastp.txt PulseNet_LT_1000_GT_250_SRRs.txt> missing_fastp_post_Rd1.txt
wc -l missing_fastp_post_Rd1.txt
#2

#Activate an environment with fastp
source /project/fsepru113/dbradshaw/miniconda3/bin/activate polish

#Run Fastp 
cat missing_fastp_post_Rd1.txt | while read read_id; do fastp --in1 reads/raw/"$read_id"_1.fastq --in2 reads/raw/"$read_id"_2.fastq --out1 reads/fastp/"$read_id".fastp_1.fastq --out2 reads/fastp/"$read_id".fastp_2.fastq --unpaired1 reads/unpaired/"$read_id".fastp.unpaired_1.fastq --unpaired2 reads/unpaired/"$read_id"_fastp.unpaired_2.fastq; done

#Check that you have right number
ls reads/fastp/*_1.fastq | wc -l 
#5763/5763

ls reads/unpaired/*_1.fastq | wc -l 
#5763/5763

ls reads/fastp/*_2.fastq | wc -l 
#5763/5763

ls reads/unpaired/*_2.fastq | wc -l 
#5763/5763


###SKESA

#Adjust step_4a_draft_skesa_assembly_reads_folder.sh to your situation
nano Step_4_PulseNet_LT_1000_GT_250_skesa_assembly.sh
sbatch Step_4_PulseNet_LT_1000_GT_250_skesa_assembly.sh


#Check for correctly made assemblies, do by size becasue will output even if it fails
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#5756/5763

#Remove assemblies that are less than 1 byte
find skesa_assembly/*.fna -type f -size -1c -exec rm "{}" \;

#Get a list of assemblies that worked the first time
ls skesa_assembly/ | cut -d "." -f 1 > assemblies_Rd1.txt
wc -l assemblies_Rd1.txt
#5756

#Make a list of unique ids of the reads that fastp worked on
grep -Fxv -f assemblies_Rd1.txt PulseNet_LT_1000_GT_250_SRRs.txt> missing_assemblies_post_Rd1.txt
wc -l missing_assemblies_post_Rd1.txt
#7

#Activate an environment with skesa
source /project/fsepru113/dbradshaw/miniconda3/bin/activate skesa

#Run Fastp 
cat missing_assemblies_post_Rd1.txt | while read read_id; do skesa --reads reads/fastp/"$read_id".fastp_1.fastq,reads/fastp/"$read_id".fastp_2.fastq --cores 4 --memory 48 > skesa_assembly/"$read_id".skesa.fna; done

#Check that you have right number
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#5763/5763

#Create a final txt to hold all SRRs that made it through the process
cp PulseNet_LT_1000_GT_250_SRRs.txt PulseNet_LT_1000_GT_250_SRRs_final.txt

#Move the files to its appropriate Set
mv *.txt *.sh Set_4_LT_1000_GT_250/

#Move the assemblies from temp folder to final folder
mv skesa_assembly/* final_assemblies
find final_assemblies/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#26696

#Remove the sras and reads but keep foldes
rm sras/*
rm reads/raw/*
rm reads/fastp/*
rm reads/unpaired/*

###Set 5 LT_250 - 5039

##Step 1 SRRs to SRAs

#Move over the txt file to main directory
mv Set_5_LT_250/PulseNet_LT_250_SRRs.txt .

#Create a srun
srun -N 1 --mem=48gb -c 24 -p short --pty $SHELL

#Get all the reads using step_1_draft_srr_to_sras.sh n= 1509
nano  Step_1_PulseNet_LT_250_SRRs.sh
sbatch  Step_1_PulseNet_LT_250_SRRs.sh

#Remove empty directories
rmdir SRR*

#Make a list of the sra files that got dounloaded
ls sras/ | cut -d "." -f 1 > Rd1_sras.txt
wc -l Rd1_sras.txt
#4999/5039

#Make a missing summarized sras file and record the number
grep -Fxv -f Rd1_sras.txt PulseNet_LT_250_SRRs.txt> missing_sras_post_downloading_Rd1.txt
wc -l missing_sras_post_downloading_Rd1.txt
#40

#Load the sra toolkit
module load sratoolkit

#Grab the missing files "locally" to see if any erros occur
prefetch --option-file missing_sras_post_downloading_Rd1.txt --verbose

#Move sras to the sra folder
cat missing_sras_post_downloading_Rd1.txt | while read srr; do mv /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/Set_5_LT_250/"$srr"/* sras/; done
rmdir SRR*

ls sras/ | wc -l
#5037/5039

#SRR169363
#SRR9660245
#No data for this SRR
#5037 is the new expected total for LT_250

#Make a new upated file based upon sras directory
ls sras/ | cut -d "." -f 1 > PulseNet_LT_250_SRRs_updated.txt
wc -l PulseNet_LT_250_SRRs_updated.txt
#5037

#Change all sras to paired reads using step_2_draft_sra_to_reads.sh
nano Step_2_PulseNet_LT_250_sras_to_reads.sh
sbatch Step_2_PulseNet_LT_250_sras_to_reads.sh

#Make a list of unique ids of the sras that got split into reads (and double check you got right number)
ls reads/raw/*_1.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R1_reads.txt
ls reads/raw/*_2.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R2_reads.txt

wc -l Rd1_R1_reads.txt
#5037

wc -l Rd1_R2_reads.txt
#5037

#No missing reads

#Run fastp on all reads without assemblies using step_3_draft_fastp_raw_reads.sh
nano Step_3_PulseNet_LT_250_fastp_raw_reads.sh
sbatch Step_3_PulseNet_LT_250_fastp_raw_reads.sh

#Make a list of unique ids of the reads that fastp worked on
ls reads/fastp/*  | cut -d "/" -f 3 | cut -d "." -f 1 | uniq > Rd1_fastp.txt
wc -l Rd1_fastp.txt
#5035/5037

#Make a missing fastp reads file and record the number
grep -Fxv -f Rd1_fastp.txt PulseNet_LT_250_SRRs_updated.txt> missing_fastp_post_Rd1.txt
wc -l missing_fastp_post_Rd1.txt
#2

#Activate an environment with fastp
source /project/fsepru113/dbradshaw/miniconda3/bin/activate polish

#Run Fastp 
cat missing_fastp_post_Rd1.txt | while read read_id; do fastp --in1 reads/raw/"$read_id"_1.fastq --in2 reads/raw/"$read_id"_2.fastq --out1 reads/fastp/"$read_id".fastp_1.fastq --out2 reads/fastp/"$read_id".fastp_2.fastq --unpaired1 reads/unpaired/"$read_id".fastp.unpaired_1.fastq --unpaired2 reads/unpaired/"$read_id"_fastp.unpaired_2.fastq; done

#Check that you have right number
ls reads/fastp/*_1.fastq | wc -l 
#5037/5037

ls reads/unpaired/*_1.fastq | wc -l 
#5037/5037

ls reads/fastp/*_2.fastq | wc -l 
#5037/5037

ls reads/unpaired/*_2.fastq | wc -l 
#5037/5037


###SKESA

#Adjust step_4a_draft_skesa_assembly_reads_folder.sh to your situation
nano Step_4_PulseNet_LT_250_skesa_assembly.sh
sbatch Step_4_PulseNet_LT_250_skesa_assembly.sh


#Check for correctly made assemblies, do by size becasue will output even if it fails
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#5027/5037

#Remove assemblies that are less than 1 byte
find skesa_assembly/*.fna -type f -size -1c -exec rm "{}" \;

#Get a list of assemblies that worked the first time
ls skesa_assembly/ | cut -d "." -f 1 > assemblies_Rd1.txt
wc -l assemblies_Rd1.txt
#5027

#Make a list of unique ids of the reads that fastp worked on
grep -Fxv -f assemblies_Rd1.txt PulseNet_LT_250_SRRs_updated.txt> missing_assemblies_post_Rd1.txt
wc -l missing_assemblies_post_Rd1.txt
#10

#Activate an environment with skesa
source /project/fsepru113/dbradshaw/miniconda3/bin/activate skesa

#Run Fastp 
cat missing_assemblies_post_Rd1.txt | while read read_id; do skesa --reads reads/fastp/"$read_id".fastp_1.fastq,reads/fastp/"$read_id".fastp_2.fastq --cores 4 --memory 48 > skesa_assembly/"$read_id".skesa.fna; done

#Check that you have right number
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#5036/5037


#Test isolate that had post-fastp reads with different number of mates
prefetch SRR10489368 --verbose
fasterq-dump -S SRR10489368/SRR10489368.sra  --verbose -e 6
fastp --in1 SRR10489368_1.fastq --in2 SRR10489368_2.fastq --out1 SRR10489368.fastp_1.fastq --out2 SRR10489368.fastp_2.fastq --unpaired1 SRR10489368.fastp.unpaired_1.fastq --unpaired2 SRR10489368_fastp.unpaired_2.fastq
#WARNNIG: different read numbers of the 2408 pack
#Read1 pack size: 106
#Read2 pack size: 1000
skesa --reads SRR10489368.fastp_1.fastq,SRR10489368.fastp_2.fastq --cores 4 --memory 48 > SRR10489368.skesa.fna
#Files SRR10489368.fastp_1.fastq,SRR10489368.fastp_2.fastq contain different number of mates

#New expected total is 5036

#Create a final txt to hold all SRRs that made it through the process
ls skesa_assembly/ | cut -d "." -f 1 > PulseNet_LT_250_SRRs_final.txt
wc -l PulseNet_LT_250_SRRs_final.txt
#5036 

#Move the files to its appropriate Set
mv *.txt *.sh Set_5_LT_250/

#Move the assemblies from temp folder to final folder
mv skesa_assembly/* final_assemblies
find final_assemblies/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#31,732

#Remove the sras and reads but keep foldes
rm sras/*
rm reads/raw/*
rm reads/fastp/*
rm reads/unpaired/*


###Senftenberg-PulseNet Project PulseNet Isolates

#Some isolates already have been assemblied in another Project (n=996)

#Get list of isolates and copy over which ones have been assemblied
cat PulseNet_PulseNet_19_23_Human_Turkey_SRRs.txt | while read assembly_id; do cp /90daydata/fsepru113/dbradshaw/Senftenberg_Projects/PulseNet_021524/all_assemblies/"$assembly_id"* ./skesa_assembly/; done

#Determine the number of isolates that were copied sucessfully
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#996/996


##Step 1 SRRs to SRAs

#Expect to have 12 isolates that had an Isolates Browser assembly
#Need to go through full process to make a SKESA assembly for these isolates

#Check for correctly made assemblies, do by size becasue will output even if it fails
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#996/996

#Get a list of assemblies that worked the first time
ls skesa_assembly/ | cut -d "." -f 1 > PN_exclusive_assemblies.txt
wc -l PN_exclusive_assemblies.txt
#996

#Make a list of unique ids of the reads that fastp worked on
grep -Fxv -f PN_exclusive_assemblies.txt PulseNet_PulseNet_19_23_Human_Turkey_SRRs.txt> Overlapping_IB_Assemblies_PN_SRRs.txt
wc -l Overlapping_IB_Assemblies_PN_SRRs.txt
#12


#Create a srun
srun -N 1 --mem=48gb -c 24 -p short --pty $SHELL

##Get sras using step_1_draft_srr_to_sras.sh scripts

#Load SRA Toolkit
module load sratoolkit

#Get all the sras in your list
prefetch --option-file Overlapping_IB_Assemblies_PN_SRRs.txt --verbose

#Then move all sra files to the new folder 
cat Overlapping_IB_Assemblies_PN_SRRs.txt | while read srr; do mv "$srr"/* sras/; done

#Remove the blank directories
rmdir SRR*

#Check that you have right number
ls sras/ | wc -l 
#12

##Step 2 SRAs to Paired Reads

##Change all sras to paired reads using scripts from step_2_draft_sra_to_reads.sh

#Split all sra files to paired reads
cat Overlapping_IB_Assemblies_PN_SRRs.txt | while read sra; do fasterq-dump -S sras/"$sra".sra -O reads/raw --verbose -e 2; done

#Check that you have right number
ls reads/raw/*_1.fastq | wc -l 
#12

##Step 3 Quality Control of Paired Reads with Fastp

##Run fastp on all paired reads using scripts from step_3_draft_fastp_raw_reads.sh

#Activate an environment with fastp
source /project/fsepru113/dbradshaw/miniconda3/bin/activate polish

#Run Fastp 
cat Overlapping_IB_Assemblies_PN_SRRs.txt | while read read_id; do fastp --in1 reads/raw/"$read_id"_1.fastq --in2 reads/raw/"$read_id"_2.fastq --out1 reads/fastp/"$read_id".fastp_1.fastq --out2 reads/fastp/"$read_id".fastp_2.fastq --unpaired1 reads/unpaired/"$read_id".fastp.unpaired_1.fastq --unpaired2 reads/unpaired/"$read_id"_fastp.unpaired_2.fastq; done

#Check that you have right number
ls reads/fastp/*_1.fastq | wc -l 
#12

ls reads/unpaired/*_1.fastq | wc -l 
#12

##Step 4 Assembly with SKESA

##Run SKESA on all paired reads using scripts from step_4a_draft_skesa_assembly_reads_folder.sh

#Activate the SKESA environment
source /project/fsepru113/dbradshaw/miniconda3/bin/activate skesa

#Run SKESA 
cat Overlapping_IB_Assemblies_PN_SRRs.txt | while read read_id; do skesa --reads reads/fastp/"$read_id".fastp_1.fastq,reads/fastp/"$read_id".fastp_2.fastq --cores 4 --memory 48 > skesa_assembly/"$read_id".skesa.fna; done

#Check that correct number was made
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#/1008

#Move the files to its appropriate Set
mv *.txt Set_6_Other_PulseNet/

#Move the assemblies from temp folder to final folder
mv skesa_assembly/* final_assemblies
find final_assemblies/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#32,740

#Remove the sras and reads but keep foldes
rm sras/*
rm reads/raw/*
rm reads/fastp/*
rm reads/unpaired/*

###Set 7 - Javiana 1014

##Step 1 SRRs to SRAs

#Create a srun
srun -N 1 --mem=48gb -c 24 -p short --pty $SHELL

#Get all the reads using step_1_draft_srr_to_sras.sh n= 1509
nano  Step_1_PulseNet_Javiana_SRRs.sh
sbatch  Step_1_PulseNet_Javiana_SRRs.sh

#Remove empty directories
rmdir SRR*

#Make a list of the sra files that got dounloaded
ls sras/ | cut -d "." -f 1 > Rd1_sras.txt
wc -l Rd1_sras.txt
#1004/1014

#Make a missing summarized sras file and record the number
grep -Fxv -f Rd1_sras.txt PulseNet_Javiana_SRRs.txt> missing_sras_post_downloading_Rd1.txt
wc -l missing_sras_post_downloading_Rd1.txt
#10

#Load the sra toolkit
module load sratoolkit

#Grab the missing files "locally" to see if any erros occur
prefetch --option-file missing_sras_post_downloading_Rd1.txt --verbose

#Move sras to the sra folder
cat missing_sras_post_downloading_Rd1.txt | while read srr; do mv /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/"$srr"/* sras/; done
rmdir SRR*

ls sras/ | wc -l
#1014/1014

#Change all sras to paired reads using step_2_draft_sra_to_reads.sh
nano Step_2_PulseNet_Javiana_sras_to_reads.sh
sbatch Step_2_PulseNet_Javiana_sras_to_reads.sh

#Make a list of unique ids of the sras that got split into reads (and double check you got right number)
ls reads/raw/*_1.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R1_reads.txt
ls reads/raw/*_2.fastq | cut -d "/" -f 3 | cut -d "_" -f 1 > Rd1_R2_reads.txt

wc -l Rd1_R1_reads.txt
#1014

wc -l Rd1_R2_reads.txt
#1014

#No missing reads

#Run fastp on all reads without assemblies using step_3_draft_fastp_raw_reads.sh
nano Step_3_PulseNet_Javiana_fastp_raw_reads.sh
sbatch Step_3_PulseNet_Javiana_fastp_raw_reads.sh

#Make a list of unique ids of the reads that fastp worked on
ls reads/fastp/*  | cut -d "/" -f 3 | cut -d "." -f 1 | uniq > Rd1_fastp.txt
wc -l Rd1_fastp.txt
#1012/1014

#Make a missing fastp reads file and record the number
grep -Fxv -f Rd1_fastp.txt PulseNet_Javiana_SRRs.txt> missing_fastp_post_Rd1.txt
wc -l missing_fastp_post_Rd1.txt
#2

#Activate an environment with fastp
source /project/fsepru113/dbradshaw/miniconda3/bin/activate polish

#Run Fastp 
cat missing_fastp_post_Rd1.txt | while read read_id; do fastp --in1 reads/raw/"$read_id"_1.fastq --in2 reads/raw/"$read_id"_2.fastq --out1 reads/fastp/"$read_id".fastp_1.fastq --out2 reads/fastp/"$read_id".fastp_2.fastq --unpaired1 reads/unpaired/"$read_id".fastp.unpaired_1.fastq --unpaired2 reads/unpaired/"$read_id"_fastp.unpaired_2.fastq; done

#Check that you have right number
ls reads/fastp/*_1.fastq | wc -l 
#1014/1014

ls reads/unpaired/*_1.fastq | wc -l 
#1014/1014

ls reads/fastp/*_2.fastq | wc -l 
#1014/1014

ls reads/unpaired/*_2.fastq | wc -l 
#1014/1014


###SKESA

#Adjust step_4a_draft_skesa_assembly_reads_folder.sh to your situation
nano Step_4_PulseNet_Javiana_skesa_assembly.sh
sbatch Step_4_PulseNet_Javiana_skesa_assembly.sh


#Check for correctly made assemblies, do by size becasue will output even if it fails
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#1012/1014

#Remove assemblies that are less than 1 byte
find skesa_assembly/*.fna -type f -size -1c -exec rm "{}" \;

#Get a list of assemblies that worked the first time
ls skesa_assembly/ | cut -d "." -f 1 > assemblies_Rd1.txt
wc -l assemblies_Rd1.txt
#1012

#Make a list of unique ids of the reads that fastp worked on
grep -Fxv -f assemblies_Rd1.txt PulseNet_Javiana_SRRs.txt> missing_assemblies_post_Rd1.txt
wc -l missing_assemblies_post_Rd1.txt
#2

#Activate an environment with skesa
source /project/fsepru113/dbradshaw/miniconda3/bin/activate skesa

#Run Fastp 
cat missing_assemblies_post_Rd1.txt | while read read_id; do skesa --reads reads/fastp/"$read_id".fastp_1.fastq,reads/fastp/"$read_id".fastp_2.fastq --cores 4 --memory 48 > skesa_assembly/"$read_id".skesa.fna; done

#Check that you have right number
find skesa_assembly/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#1014/1014

#Create a final txt to hold all SRRs that made it through the process
ls skesa_assembly/ | cut -d "." -f 1 > PulseNet_Javiana_SRRs_final.txt
wc -l PulseNet_Javiana_SRRs_final.txt
#1014 

#Move the files to its appropriate Set
mv *.txt *.sh Set_7_Javiana/

#Move the assemblies from temp folder to final folder
mv skesa_assembly/* final_assemblies
find final_assemblies/*.fna -type f -size +1c -exec ls "{}" \; | wc -l
#33,754

#Remove the sras and reads but keep foldes
rm sras/*
rm reads/raw/*
rm reads/fastp/*
rm reads/unpaired/*

#Make a list of all assemblies
ls final_assemblies/ | cut -d "." -f 1 > PulseNet_final_assemblies.txt
wc -l PulseNet_final_assemblies.txt
#33754

##Step 5 Bakta Annotation

#Make a new directory
mkdir bakta_full_annotation

#Adjust step_5_draft_bakta_annotation.sh for first 10k isolates
nano Step_5_PulseNet_All_Bakta_Annotation_Pt1.sh
sbatch Step_5_PulseNet__All_Bakta_Annotation_Pt1.sh

#If above 10000, then use Lines 14-15 of step_5_draft_bakta_full_annotation.sh for first 10000
#Then adjust comment out those lines, and edit and use Lines 18-19 for the rest
#Keep in mind to adjust the count accordingly in line 18 for each new 10000
#Line 4 -> #SBATCH --array=1-4062
#Line 18 -> proteome_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_bakta_protein_faas/*.faa | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#Adjustments above are acually running proteomes 10001-14062
nano Step_5_PulseNet_All_Bakta_Annotation_Pt2.sh
sbatch Step_5_PulseNet_All_Bakta_Annotation_Pt2.sh

#Double check correct number was made
ls bakta_full_annotation/*/*.txt | wc -l 

#Make a list of assemblies that were annotated correctly
cd bakta_full_annotation
ls */*.txt | cut -d "/" -f 2 | cut -d "." -f 1 > ../full_bakta_annotation_post_Pt1_Rd1.txt
wc ../full_bakta_annotation_post_Pt1_Rd1.txt -l 
#33733

#Make a list of missing full anotation files files
cd ..
grep -Fxv -f full_bakta_annotation_post_Pt1_Rd1.txt PulseNet_final_assemblies.txt> missing_full_bakta_anno_Pt1_Rd1.txt
wc missing_full_bakta_anno_Pt1_Rd1.txt -l
#21

#Run locally to see any errors

#Activate the environment
source /project/fsepru113/dbradshaw/miniconda3/bin/activate bakta

cat missing_full_bakta_anno_Pt1_Rd1.txt | while read sample_id; do bakta --db bakta --db /project/fsepru113/dbradshaw/databases/bakta_db_full/ -o bakta_full_annotation/$sample_id.bakta --keep-contig-headers -t 24 --genus Salmonella --species enterica --verbose --force --prefix $sample_id /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/final_assemblies/$sample_id*; done

#Doublecheck they populated correctly
ls bakta_full_annotation/*/*.txt | wc -l 
#33754

#Make a directory to copy the protein files over to
mkdir /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_bakta_protein_faas

#Move bakta faa files to appropriate folder by editing and running draft_mv_or_cp_array.sh Pts 1-x
nano Bakta_faas_to_Epitope_Blasting_Pt1.sh
sbatch Bakta_faas_to_Epitope_Blasting_Pt1.sh

#Adjust draft_mv_or_cp_array.sh for next 10k isolates and repeat as necessary
nano Bakta_faas_to_Epitope_Blasting_Pt2.sh
sbatch Bakta_faas_to_Epitope_Blasting_Pt2.sh

#Check you have correct number of faa files
#May have to do following steps in the directory itself due to this error - -bash: /usr/bin/wc: Argument list too long, adjust accordingly
wc /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_bakta_protein_faas/* -l 
#33753

#Make a list of those that were copied over
ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_bakta_protein_faas/ | cut -d "." -f 1 > mv_bakta_faas_Rd1.txt

#Make a list of ones that were not copied over
grep -Fxv -f mv_bakta_faas_Rd1.txt PulseNet_final_assemblies.txt> missing_mv_bakta_faas_Rd1.txt
cat missing_mv_bakta_faas_Rd1.txt
#SRR9999933

#cd into that directory, make sure bakta worked, and then try to recopy it over
cd bakta_full_annotation/SRR9999933.bakta/
cp SRR9999933.faa /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_bakta_protein_faas/

#Check you have correct number of faa files
find /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_bakta_protein_faas/* -type f -size +1c -exec ls "{}" \; | wc -l
wc /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_bakta_protein_faas/* -l 
#33754


##CheckM QC

#Adjust step_qc_draft_checkm.sh as needed
nano Step_QC_PulseNet_Set_1_checkM.sh
sbatch Step_QC_PulseNet_Set_1_checkM.sh 

#Ran as Sets to ensure quicker analysis

#Concatenate all bin_stats_ext.tsvs
cat checkm/*/storage/bin_stats_ext.tsv > PulseNet_checkm_results.tsv

#Send down to Excel and edit
#First select the first column and Use Text to Column to separate by "," and ":"
#Then delete everthing before "Completeness" (B-U)
#Delete everything after and including   'Translation table' (AF-BLW)
#Delete the non-numerical columns i.e. column names
#Add a row and copy in the following column names, should separate automatically due to commas
#Genome.name,Completeness,Contamination,GC,GC.std,Genome.size,Ambiguous.bases,Scaffolds,Contigs,Longest.scaffold,Longest.contig,N50.scaffolds,N50.contigs,Mean.scaffold.length,Mean.contig.length,Coding.density
#If needed remove the ".skesa" from the Genome.name column values

###Epitope Blasting vs Proteomes - Steps 6 and 7###


#Make the following directories for the epitope blasting steps
cd /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/

mkdir V1
mv * V1/
mv V1/Pre_Positive_Homology_UK1_Epitopes.fasta .
mv V1/blast_dbs/ .
mv V1/protein_data/ .

mkdir blast_dbs
mkdir epitope_blasting_results
mkdir stds_epitope_blasting_results
mkdir epitope_blasting_results_summaries_Step_1
mkdir stds_epitope_blasting_results_summaries_Step_1
mkdir PulseNet_splitting
mkdir epitope_blasting_results_summaries_Step_2
mkdir stds_epitope_blasting_results_summaries_Step_2
mkdir epitope_blasting_results_summaries_Step_3
mkdir stds_epitope_blasting_results_summaries_Step_3
mkdir epitope_blasting_results_summaries_Step_4
mkdir stds_epitope_blasting_results_summaries_Step_4


#Copy over the .Renviron file to that R can work
 cp ./project/fsepru113/dbradshaw/.Renviron .

#Copy the fasta file of the construct epitopes to HPC
nano VaxiJen_Localization_Agnostic_Construct_Epitopes.fasta

##Step 6 Proteome Blasting

#Adjust step_6_draft_positive_homology_epitope_proteome_blasting.sh for first 10k isolates
nano Step_6_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Pt1.sh
sbatch Step_6_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Pt1.sh

#Adjust step_6_draft_positive_homology_epitope_proteome_blasting.sh for next 10k isolates and repeat as necessary
nano Step_6_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Pt2.sh
sbatch Step_6_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Pt2.sh

#Double check correct number was made
ls blast_dbs/*.pos | wc -l 
ls epitope_blasting_results | wc -l

#Make a list of assemblies that were blasted correctly
ls epitope_blasting_results/*.txt | cut -d "/" -f 2 | cut -d "_" -f 1 > proteome_blasting_post_Rd1.txt
wc proteome_blasting_post_Rd1.txt -l

#Make a list of missing full anotation files files
grep -Fxv -f proteome_blasting_post_Rd1.txt PulseNet_final_assemblies.txt > missing_proteome_blasting_Rd1.txt
wc missing_proteome_blasting_Rd1 -l
#GCA_020652735.1

#Run with step_6_draft_positive_homology_epitope_proteome_blasting_redo.sh
nano Step_6_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Redo_Pt1.sh
sbatch Step_6_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Redo_Pt1.sh

#Doublecheck they populated correctly
ls blast_dbs | wc -l 
ls epitope_blasting_results | wc -l

#Edit the Step_7a_draft_proteome_blasting_summaries_Step_1_R.sh
nano PulseNet_proteome_blasting_summaries_Step_1.R
 
#Edit the step_7a_draft_positive_homology_epitope_proteome_blasting_summaries_Step_1.sh as needed
nano Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Pt1.sh
sbatch Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Pt1.sh

#Adjust step_7a_draft_positive_homology_epitope_proteome_blasting_summaries_Step_1.sh for next 10k isolates and repeat as necessary
nano Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Pt2.sh
sbatch Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Pt2.sh

#Determine what proteomes were actually summarized and record the number of proteomes
find epitope_blasting_results_summaries_Step_1/ -type f -size +1c -exec ls "{}" \; | wc -l
#PulseNet - 33754/33754

#If you do not have all the summaries then do the following steps, if you do then skip to "cd epitope_blasting_results_summaries_Step_1/"

#Make a list of the summarized proteomes
ls epitope_blasting_results_summaries_Step_1/ | cut -d "_" -f 1 > epitope_blasting_results_summaries_Step_1_Rd1.txt

#Make a missing summarized proteomes file and record the number
grep -Fxv -f epitope_blasting_results_summaries_Step_1_Rd1.txt PulseNet_final_assemblies.txt > missing_epitope_blasting_results_summaries_Step_1_Rd1.txt

#Run with step_7a_draft_positive_homology_epitope_proteome_blasting_summaries_Step_1_redo.sh
nano Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Redo_Pt1.sh
sbatch Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Redo_Pt1.sh

#Double check you now have all the necessary summaries and record number
ls epitope_blasting_results_summaries_Step_1/ | wc -l
#PulseNet - 33754

#If you do then cd into that directory, if not then redo above steps
cd epitope_blasting_results_summaries_Step_1/

#Get the name of the first file in the folder and use in next script
ls * | head -n 1
#SRR10000867_proteome_epitope_homology.txt

#Concatenate all summaries by first grabbing header from one and then results from all
head -n 1 SRR10000867_proteome_epitope_homology.txt > PulseNet_proteomes_epitope_homology.txt; tail -n +2 -q *.txt >> PulseNet_proteomes_epitope_homology.txt

#If need be, you can also use edit and use Draft_Proteome_Blasting_Summaries_Post_Step_1_Concatenation.sh
nano PulseNet_combining_summaries.sh
sbatch PulseNet_combining_summaries.sh
#Number of lines should be about number of epitopes x number of downloaded proteomes
#Theoretical PulseNet - 28 x 33754 = 945,112

wc -l PulseNet_proteomes_epitope_homology.txt
#Actual PulseNet - 945,045

#Move it to main directory 
mv PulseNet_proteomes_epitope_homology.txt ../
cd ../

#Copy down to local environment to remove isolates that did not pass the CheckM filter
scp david.bradshaw@ceres.scinet.usda.gov:/project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_proteomes_epitope_homology.txt .

#Copy backup to HPC for analysis
scp PulseNet_proteomes_epitope_homology_Step_1_QCed.txt david.bradshaw@ceres.scinet.usda.gov:/project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/

#Copy to the PulseNet_splitting directory
cp PulseNet_proteomes_epitope_homology_Step_1_QCed.txt PulseNet_splitting/

#Change into this directory
cd PulseNet_splitting

#Split PulseNet_proteomes_epitope_homology_Step_1_QCed.txt into a summary of each epitope (Unique_Id column = column 1)
awk '{print > $1".txt"}' PulseNet_proteomes_epitope_homology_Step_1_QCed.txt

#Remove the extra files
rm PulseNet_proteomes_epitope_homology_Step_1_QCed.txt
rm Unique_Id.txt

#Check to make sure there is the correct number of files i.e. matches number of non-unique epitopes
ls | wc -l #28

#Adjust Step_7b_draft_proteome_blasting_summaries_Step_2_R.sh and step_7b_draft_positive_homology_epitope_proteome_blasting_summaries_Step_2.sh as needed
nano Step_7b_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_2_Pt1.sh
sbatch Step_7b_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_2_Pt1.sh 
#Informs Supplementary Table 12

#Determine what epitopes were actually summarized and record the number of epitopes
ls epitope_blasting_results_summaries_Step_2/ | wc -l
#PulseNet - 28

#Move into the directory and combine all files into one
cd epitope_blasting_results_summaries_Step_2

ls | head -n 1
#AEF05957.1_MHCII_Pos_337_proteomes_epitope_homology_summary.txt


head -n 1 AEF05957.1_MHCII_Pos_337_proteomes_epitope_homology_summary.txt > PulseNet_proteomes_epitope_homology_summary.tsv; tail -n +2 -q *.txt >> PulseNet_proteomes_epitope_homology_summary.tsv

#Move it to the main directory
mv PulseNet_proteomes_epitope_homology_summary.tsv ../PulseNet_proteomes_epitope_homology_summary.txt

#Copy down to local environment to remove isolates that did not pass the CheckM filter
scp david.bradshaw@ceres.scinet.usda.gov:/project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_proteomes_epitope_homology_summary.txt .



#Adjust Step_7c_draft_proteome_blasting_summaries_Step_3_R.sh and step_7c_draft_positive_homology_epitope_proteome_blasting_summaries_Step_3.sh as needed
nano Step_7c_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_3_Pt1.sh
sbatch Step_7c_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_3_Pt1.sh 

#Determine what epitopes were actually summarized and record the number of epitopes
ls epitope_blasting_results_summaries_Step_3/ | wc -l
#PulseNet - 28

#Move into the directory and combine all files into one
cd epitope_blasting_results_summaries_Step_3

ls | head -n 1
#AEF05957.1_MHCII_Pos_337_proteomes_epitope_homology_summary.txt


head -n 1 AEF05957.1_MHCII_Pos_337_proteomes_epitope_homology_summary.txt > PulseNet_proteomes_epitope_homology_summary.tsv; tail -n +2 -q *.txt >> PulseNet_proteomes_epitope_homology_summary.tsv

awk 'NR==1 { print; next } FNR > 1' *.tsv > PulseNet_proteomes_epitope_homology_serovar_summary.tsv


#Move it to the main directory
mv PulseNet_proteomes_epitope_homology_summary.tsv ../PulseNet_proteomes_epitope_homology_summary.txt

#Copy down to local environment to remove isolates that did not pass the CheckM filter
scp david.bradshaw@ceres.scinet.usda.gov:/project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_proteomes_epitope_homology_summary.txt .

#Adjust Step_7d_draft_proteome_blasting_summaries_Step_4_R.sh and step_7d_draft_positive_homology_epitope_proteome_blasting_summaries_Step_4.sh as needed
nano Step_7d_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_4_Pt1.sh
sbatch Step_7d_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_4_Pt1.sh 
#Informs Supplementary Table 12

#Determine what epitopes were actually summarized and record the number of epitopes
ls epitope_blasting_results_summaries_Step_4/ | wc -l
#PulseNet - 28

#Move into the directory and combine all files into one
cd epitope_blasting_results_summaries_Step_4

awk 'NR==1 { print; next } FNR > 1' *.tsv > ../PulseNet_proteomes_epitope_homology_serovar_counts.tsv

#Copy down to local environment to remove isolates that did not pass the CheckM filter
scp david.bradshaw@ceres.scinet.usda.gov:/project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_proteomes_epitope_homology_summary.txt .

#Adjust Step_7e_draft_proteome_blasting_summaries_Step_5_R.sh and step_7e_draft_positive_homology_epitope_proteome_blasting_summaries_Step_5.sh as needed
nano Step_7e_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_5_Pt1.sh
sbatch Step_7e_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_5_Pt1.sh 

#Determine what epitopes were actually summarized and record the number of epitopes
ls epitope_blasting_results_summaries_Step_5/ | wc -l
#PulseNet - 28

#Move into the directory and combine all files into one
cd epitope_blasting_results_summaries_Step_5

awk 'NR==1 { print; next } FNR > 1' *.tsv > ../PulseNet_proteomes_epitope_homology_isosrc_summary.tsv

#Adjust Step_7f_draft_proteome_blasting_summaries_Step_6_R.sh and step_7f_draft_positive_homology_epitope_proteome_blasting_summaries_Step_6.sh as needed
nano Step_7f_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_6_Pt1.sh
sbatch Step_7f_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_6_Pt1.sh 
#Informs Supplementary Table 13

#Determine what epitopes were actually summarized and record the number of epitopes
ls epitope_blasting_results_summaries_Step_6/ | wc -l
#PulseNet - 28

#Move into the directory and combine all files into one
cd epitope_blasting_results_summaries_Step_6

awk 'NR==1 { print; next } FNR > 1' *.tsv > ../PulseNet_proteomes_epitope_homology_isosrc_counts.tsv

#Adjust Step_7g_draft_proteome_blasting_summaries_Step_7_R.sh and step_7g_draft_positive_homology_epitope_proteome_blasting_summaries_Step_7.sh as needed
nano Step_7g_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_7_Pt1.sh
sbatch Step_7g_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_7_Pt1.sh 

#Determine what epitopes were actually summarized and record the number of epitopes
ls epitope_blasting_results_summaries_Step_7/ | wc -l
#PulseNet - 28

#Move into the directory and combine all files into one
cd epitope_blasting_results_summaries_Step_7

awk 'NR==1 { print; next } FNR > 1' *.tsv > ../PulseNet_proteomes_epitope_homology_otbksrc_summary.tsv

#Adjust Step_7h_draft_proteome_blasting_summaries_Step_8_R.sh and step_7h_draft_positive_homology_epitope_proteome_blasting_summaries_Step_8.sh as needed
nano Step_7h_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_8_Pt1.sh
sbatch Step_7h_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_8_Pt1.sh 
#Informs Supplementary Table 14

#Determine what epitopes were actually summarized and record the number of epitopes
ls epitope_blasting_results_summaries_Step_8/ | wc -l
#PulseNet - 28

#Move into the directory and combine all files into one
cd epitope_blasting_results_summaries_Step_8

awk 'NR==1 { print; next } FNR > 1' *.tsv > ../PulseNet_proteomes_epitope_homology_otbksrc_counts.tsv

#Copy the following files from SCINet down to you local computer into the appropriate folder by changing directories on your local computer to target location and doing scp
#PulseNet_proteomes_epitope_homology_summary.txt
#PulseNet_proteomes_epitope_homology_Step_1_QCed.txt
#Delete what you believe you do not need (the stdout and stderr likely won't have much info for instance)

#Head over to R for summarization of the number of epitopes with 100% identity with varying percentages of proteomes per serovar

#...

#Importing information from R & Excel...
###69 Non-Unique MHCI with 100% Identity to 99% of All Serovars' Proteomes###
###216 Non-Unique MHCII with 100% Identity to 99% of All Serovars' Proteomes###

#...
















