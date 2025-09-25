#!/bin/bash
#SBATCH --job-name=Enteritidis_proteome_blasting_summaries
#SBATCH -p short
#SBATCH --array=1-15
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"

blast_file="$(ls Enteritidis_missing_summaries/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
proteome_name="$(echo $blast_file | cut -d "_" -f 1,2)"

#Load the module
module load r/4.3.0

#Run the RScript
Rscript Enteritidis_proteome_blasting_summaries_Step_1.R $blast_file $proteome_name