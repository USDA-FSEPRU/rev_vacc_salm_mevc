#!/bin/bash
#SBATCH --job-name=Step_7g_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_7
#SBATCH -p ceres
#SBATCH --array=3-28
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Step_7g_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_7.stdout.%j.%N"
#SBATCH -e "Step_7g_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_7.stderr.%j.%N" 
#SBATCH --account=fsepru113

#Create a variables to store file and epitope names (first 10k set)
summary_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_splitting/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
epitope_unique_id="$(echo $summary_file | cut -d "." -f 1,2)"

#Load the module
module load r/4.3.0

#Run the RScript
Rscript PulseNet_proteome_blasting_summaries_Step_7.R $summary_file $epitope_unique_id
