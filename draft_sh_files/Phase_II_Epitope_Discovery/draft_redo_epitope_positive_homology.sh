#!/bin/bash
#SBATCH --job-name=Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt1_V2
#SBATCH -p short
#SBATCH --array=1-10000
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stdout.%j.%N"
#SBATCH -e "Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Enteritidis_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Enteritidis_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name

#!/bin/bash
#SBATCH --job-name=Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt2_V2
#SBATCH -p short
#SBATCH --array=1-10000
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt2_V2.stdout.%j.%N"
#SBATCH -e "Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt2_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Enteritidis_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Enteritidis_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name

#!/bin/bash
#SBATCH --job-name=Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt3_V2
#SBATCH -p short
#SBATCH --array=1-10000
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt3_V2.stdout.%j.%N"
#SBATCH -e "Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt3_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Enteritidis_faas/epitope_blasting_results/ | head -n 299115 | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Enteritidis_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+20000)) | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name

#!/bin/bash
#SBATCH --job-name=Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt4_V2
#SBATCH -p short
#SBATCH --array=1-7000
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt4_V2.stdout.%j.%N"
#SBATCH -e "Enteritidis_Proteomes_Blasting_Summaries_Step_1_Pt4_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Enteritidis_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Enteritidis_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+30000)) | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name


#!/bin/bash
#SBATCH --job-name=Typhimurium_Proteomes_Blasting_Summaries_Step_1_Pt1_V2
#SBATCH -p short
#SBATCH --array=1-10000
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Typhimurium_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stdout.%j.%N"
#SBATCH -e "Typhimurium_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Typhimurium_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Typhimurium_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name

#!/bin/bash
#SBATCH --job-name=Typhimurium_Proteomes_Blasting_Summaries_Step_1_Pt2_V2
#SBATCH -p short
#SBATCH --array=1-10000
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Typhimurium_Proteomes_Blasting_Summaries_Step_1_Pt2_V2.stdout.%j.%N"
#SBATCH -e "Typhimurium_Proteomes_Blasting_Summaries_Step_1_Pt2_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Typhimurium_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Typhimurium_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name

#!/bin/bash
#SBATCH --job-name=Typhimurium_Proteomes_Blasting_Summaries_Step_1_Pt3_V2
#SBATCH -p short
#SBATCH --array=1-5983
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Typhimurium_Proteomes_Blasting_Summaries_Step_1_Pt3_V2.stdout.%j.%N"
#SBATCH -e "Typhimurium_Proteomes_Blasting_Summaries_Step_1_Pt3_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Typhimurium_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Typhimurium_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+20000)) | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name











#!/bin/bash
#SBATCH --job-name=Infantis_Proteomes_Blasting_Summaries_Step_1_Pt1_V2
#SBATCH -p short
#SBATCH --array=1-10000
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Infantis_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stdout.%j.%N"
#SBATCH -e "Infantis_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Infantis_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Infantis_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name

#!/bin/bash
#SBATCH --job-name=Infantis_Proteomes_Blasting_Summaries_Step_1_Pt2_V2
#SBATCH -p short
#SBATCH --array=1-4052
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Infantis_Proteomes_Blasting_Summaries_Step_1_Pt2_V2.stdout.%j.%N"
#SBATCH -e "Infantis_Proteomes_Blasting_Summaries_Step_1_Pt2_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Infantis_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Infantis_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name













#!/bin/bash
#SBATCH --job-name=Kentucky_Proteomes_Blasting_Summaries_Step_1_Pt1_V2
#SBATCH -p short
#SBATCH --array=1-10000
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Kentucky_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stdout.%j.%N"
#SBATCH -e "Kentucky_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Kentucky_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Kentucky_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name

#!/bin/bash
#SBATCH --job-name=Kentucky_Proteomes_Blasting_Summaries_Step_1_Pt2_V2
#SBATCH -p short
#SBATCH --array=1-1100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Kentucky_Proteomes_Blasting_Summaries_Step_1_Pt2_V2.stdout.%j.%N"
#SBATCH -e "Kentucky_Proteomes_Blasting_Summaries_Step_1_Pt2_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Kentucky_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Kentucky_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name







#!/bin/bash
#SBATCH --job-name=Hadar_Proteomes_Blasting_Summaries_Step_1_Pt1_V2
#SBATCH -p short
#SBATCH --array=1-1704
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Hadar_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stdout.%j.%N"
#SBATCH -e "Hadar_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Hadar_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Hadar_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name











#!/bin/bash
#SBATCH --job-name=Uganda_Proteomes_Blasting_Summaries_Step_1_Pt1_V2
#SBATCH -p short
#SBATCH --array=1-961
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Uganda_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stdout.%j.%N"
#SBATCH -e "Uganda_Proteomes_Blasting_Summaries_Step_1_Pt1_V2.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Uganda_faas/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Uganda_faas/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1,2)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript Proteome_blasting_summaries_Step_1_V2.R $blast_file $proteome_name





ls epitope_blasting_results_summaries_Step_1/ | cut -d "_" -f 1,2 > epitope_blasting_results_summaries_Step_1_Rd1.txt

#Make a missing summarized proteomes file and record the number
grep -Fxv -f epitope_blasting_results_summaries_Step_1_Rd1.txt Kentucky_assemblies.txt > missing_epitope_blasting_results_summaries_Step_1_Rd1.txt

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

head -n 1 GCA_000009505.1_proteome_epitope_homology.txt > Enteritidis_proteomes_epitope_homology.txt; tail -n +2 -q *.txt >> Enteritidis_proteomes_epitope_homology.txt
head -n 1 GCA_000006945.1_proteome_epitope_homology.txt > Typhimurium_proteomes_epitope_homology.txt; tail -n +2 -q *.txt >> Typhimurium_proteomes_epitope_homology.txt
head -n 1 GCA_000230875.1_proteome_epitope_homology.txt > Infantis_proteomes_epitope_homology.txt; tail -n +2 -q *.txt >> Infantis_proteomes_epitope_homology.txt
head -n 1  > Kentucky_proteomes_epitope_homology.txt; tail -n +2 -q *.txt >> Kentucky_proteomes_epitope_homology.txt
head -n 1 GCA_000171515.1_proteome_epitope_homology.txt > Hadar_proteomes_epitope_homology.txt; tail -n +2 -q *.txt >> Hadar_proteomes_epitope_homology.txt
head -n 1 GCA_000231705.2_proteome_epitope_homology.txt > Uganda_proteomes_epitope_homology.txt; tail -n +2 -q *.txt >> Uganda_proteomes_epitope_homology.txt

#Number of lines should be about number of epitopes x number of downloaded proteomes
#Theoretical Enteritidis - 526 x 37000 = 19462000
#Theoretical Hadar - 526 x 1704 = 896304
#Theoretical Infantis - 526 x 14052 = 7391352
#Theoretical Kentucky - 526 x 11100 = 5838600
#Theoretocal Typhimurium - 526 x 25983 = 13667058
#Theoretical Uganda - 526 x 961 = 505486


wc -l PulseNet_proteomes_epitope_homology.txt
#Actual Enteritidis - 19,461,856
#Actual Typhimurium - 945,045
#Actual Infantis - 945,045
#Actual Kentucky - 945,045
#Actual Hadar - 945,045
#Actual Uganda - 945,045


#Move it to main directory and copy to the PulseNet_splitting directory
mv Enteritidis_proteomes_epitope_homology.txt ../
cd ../
cp Enteritidis_proteomes_epitope_homology.txt Enteritidis_splitting/

Change into this directory
cd Enteritidis_splitting

#Split Enteritidis_proteomes_epitope_homology.txt into a summary of each epitope (Unique_Id column = column 1)
awk '{print > $1".txt"}' Enteritidis_proteomes_epitope_homology.txt

#Remove the extra files
rm Enteritidis_proteomes_epitope_homology.txt
rm Unique_Id.txt

#Check to make sure there is the correct number of files i.e. matches number of non-unique epitopes
ls | wc -l #526





#!/bin/bash
#SBATCH --job-name=Kentucky_Proteomes_Blasting_Summaries_Step_2_V2
#SBATCH -p short
#SBATCH --array=1-526
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Kentucky_Proteomes_Blasting_Summaries_Step_2_V2.stdout.%j.%N"
#SBATCH -e "Kentucky_Proteomes_Blasting_Summaries_Step_2_V2.stderr.%j.%N" 

#Create a variables to store file and epitope names (first 10k set)
summary_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Kentucky_faas/Kentucky_splitting/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
epitope_unique_id="$(echo $summary_file | cut -d "." -f 1,2)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#summary_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/UK1_Epitope_Proteome_Blasting/Kentucky_faas/Kentucky_splitting/ | head -n ((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#epitope_unique_id="$(echo $summary_file | cut -d "." -f 1,2)"

#Load the module
module load r/4.3.0

#Run the RScript
Rscript Kentucky_proteome_blasting_summaries_Step_2_V2.R $summary_file $epitope_unique_id








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
path <- paste0("Kentucky_splitting/", args[6])

#Upload results from summarizing blasting proteome results for one epitope
proteomes_epitope_homology <- read.delim(path, header=FALSE, sep=" ")

##100% Percent Identity with Varying Coverage

#Determine percentage of proteomes with 100% pident
pident_epitope_percentage <- sum(proteomes_epitope_homology$V3=="TRUE")/11100*100   #Calculate the percentage of TRUEs across all proteomes

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
qcovhsp_epitope_percentage <- sum(proteomes_epitope_homology$V5=="TRUE")/11100*100   #Calculate the percentage of TRUEs across all proteomes

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
full_0_epitope_percentage <- sum(proteomes_epitope_homology$V6=="TRUE")/11100*100   #Calculate the percentage of TRUEs across all proteomes

#Create a series of TRUE/FALSE statements regarding how high this percentage is
full_0_ET_100 <- ifelse(full_0_epitope_percentage == 100, TRUE, FALSE) #If the percentage is equal to 100, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_99 <- ifelse(full_0_epitope_percentage > 99, TRUE, FALSE) #If the percentage is greater than to 99, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_98 <- ifelse(full_0_epitope_percentage > 98, TRUE, FALSE) #If the percentage is greater than to 98, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_97 <- ifelse(full_0_epitope_percentage > 97, TRUE, FALSE) #If the percentage is greater than to 97, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_96 <- ifelse(full_0_epitope_percentage > 96, TRUE, FALSE) #If the percentage is greater than to 96, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_95 <- ifelse(full_0_epitope_percentage > 95, TRUE, FALSE) #If the percentage is greater than to 95, say homology is TRUE, otherwise say Homology is FALSE
full_0_GT_90 <- ifelse(full_0_epitope_percentage > 90, TRUE, FALSE) #If the percentage is greater than to 90, say homology is TRUE, otherwise say Homology is FALSE

#Determine percentage of proteomes with >88% coverage and percent identity
partial_1_epitope_percentage <- sum(proteomes_epitope_homology$V7=="TRUE")/11100*100   #Calculate the percentage of TRUEs across all proteomes

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
colnames(proteomes_epitope_homology_summary) <- c("Unique_Id", "Epitope_vs_Kentucky_Pident_Percentage", "Pident_Kentucky_ET_100", "Pident_Kentucky_GT_99", "Pident_Kentucky_GT_98", "Pident_Kentucky_GT_97", "Pident_Kentucky_GT_96", "Pident_Kentucky_GT_95", "Pident_Kentucky_GT_90",
                                                  "Epitope_vs_Kentucky_Qcovhsp_Percentage", "Qcovhsp_Kentucky_ET_100", "Qcovhsp_Kentucky_GT_99", "Qcovhsp_Kentucky_GT_98", "Qcovhsp_Kentucky_GT_97", "Qcovhsp_Kentucky_GT_96", "Qcovhsp_Kentucky_GT_95", "Qcovhsp_Kentucky_GT_90",
                                                  "Epitope_vs_Kentucky_Full_0_Percentage", "Full_0_Kentucky_ET_100", "Full_0_Kentucky_GT_99", "Full_0_Kentucky_GT_98", "Full_0_Kentucky_GT_97", "Full_0_Kentucky_GT_96", "Full_0_Kentucky_GT_95", "Full_0_Kentucky_GT_90",
                                                  "Epitope_vs_Kentucky_Partial_1_Percentage", "Partial_1_Kentucky_ET_100", "Partial_1_Kentucky_GT_99", "Partial_1_Kentucky_GT_98", "Partial_1_Kentucky_GT_97", "Partial_1_Kentucky_GT_96", "Partial_1_Kentucky_GT_95", "Partial_1_Kentucky_GT_90")

#Split Unique_Id column into Identity and Extra bits
proteomes_epitope_homology_summary <- separate(data = proteomes_epitope_homology_summary, col = Unique_Id, into = c("Identity", "Extra"), sep = "_", extra = "merge", remove = FALSE)

#Write out the summary for this epitope
write_delim(proteomes_epitope_homology_summary, paste0("epitope_blasting_results_summaries_Step_2/", Unique_Id, "_proteomes_epitope_homology_summary.txt"))