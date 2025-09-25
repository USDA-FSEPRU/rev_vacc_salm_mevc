#!/bin/bash
#SBATCH --job-name=Step_7b_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_2_Redo
#SBATCH -p short
#SBATCH --array=1-2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "Step_7b_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_2.stdout.%j.%N"
#SBATCH -e "Step_7b_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Summaries_Step_2.stderr.%j.%N" 

#Create a variables to store file and epitope names (first 10k set)

epitope_unique_id="$(cat missing_epitope_blasting_results_summaries_Step_2_Rd1.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
summary_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_splitting/$epitope_unique_id* | cut -d "/" -f 8)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#summary_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_splitting/ | head -n ((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#epitope_unique_id="$(echo $summary_file | cut -d "." -f 1,2)"

#Load the module
module load r/4.3.0

#Run the RScript
Rscript PulseNet_proteome_blasting_summaries_Step_2.R $summary_file $epitope_unique_id