#!/bin/bash
#SBATCH --job-name=Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Redo_Pt1
#SBATCH -p short
#SBATCH --array=1-2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Redo_Pt1.stdout.%j.%N"
#SBATCH -e "Step_7a_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_1_Redo_Pt1.stderr.%j.%N" 

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
proteome_name="$(cat missing_epitope_blasting_results_summaries_Step_1_Rd1.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/epitope_blasting_results/$proteome_name* | cut -d "/" -f 8)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#proteome_name="$(cat missing_epitope_blasting_results_summaries_Step_1_Rd1.txt | head -n $((SLURM_ARRAY_TASK_ID+10000) | tail -n 1)"
#blast_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/epitope_blasting_results/$proteome_name* | cut -d "/" -f 8)"

#Load the module
module load r/4.3.0

#Run the RScript
Rscript PulseNet_proteome_blasting_summaries_Step_1.R $blast_file $proteome_name