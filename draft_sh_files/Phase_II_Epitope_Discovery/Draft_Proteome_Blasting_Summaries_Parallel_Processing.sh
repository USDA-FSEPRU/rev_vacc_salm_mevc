#!/bin/bash
#SBATCH --job-name=Kentucky_blast_results_summarizing
#SBATCH -p medium
#SBATCH --mem-per-cpu=12gb
#SBATCH -c 24
#SBATCH -o "Kentucky_blast_results_summarizing.stdout.%j.%N"
#SBATCH -e "Kentucky_blast_results_summarizing.stderr.%j.%N"

#Load the module
module load r/4.3.0

#Run the Rscript
Rscript Kentucky_blast_results_summarizing.R