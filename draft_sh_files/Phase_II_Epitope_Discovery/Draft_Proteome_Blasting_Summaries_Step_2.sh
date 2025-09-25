#!/bin/bash
#SBATCH --job-name=Infantis_proteome_blasting_summaries
#SBATCH -p short
#SBATCH --array=1-2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 1:00:00
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N" 

summary_file="$(ls Infantis_splitting/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
epitope_unique_id="$(echo $summary_file | cut -d "." -f 1,2)"

#Load the module
module load r/4.3.0

#Run the RScript
Rscript Infantis_blast_results_summarizing_Step_2_array.R $summary_file $epitope_unique_id