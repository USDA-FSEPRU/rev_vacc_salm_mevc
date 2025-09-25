#!/bin/bash
#SBATCH --job-name=Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Pt1
#SBATCH -p short
#SBATCH --array=1-10
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Pt1.stdout.%j.%N"
#SBATCH -e "Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Pt1.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/epitope_blasting_results/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/epitope_blasting_results/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#proteome_name="$(echo $blast_file | cut -d "/" -f 8  | cut -d "_" -f 1)"


#Load the module
module load r/4.3.0

#Run the RScript
Rscript PulseNet_proteome_blasting_summaries_Step_1.R $blast_file $proteome_name