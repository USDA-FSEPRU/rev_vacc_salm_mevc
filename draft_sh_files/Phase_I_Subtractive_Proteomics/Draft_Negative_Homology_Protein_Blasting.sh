#!/bin/bash
#SBATCH --job-name=negative_protein_blast_Sus_scrofa
#SBATCH -p short
#SBATCH --mem=12gb
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -o "Sus_scrofa.protein.negative.blasting.stdout.%j.%N"
#SBATCH -e "Sus_scrofa.protein.negative.blasting.stderr.%j.%N"

source /project/fsepru113/dbradshaw/miniconda3/bin/activate blast

#Change the strain name to your different targets as need

#Blast against the remote nr database but focus on your organism of interest
#Create an archive output to convert to other formats as needed
blastp -query Post_VaxiJen_UK1_Protein_Accessions.txt -remote -db nr -outfmt 11 -out UK1_Sus_scrofa_archive_070323.txt -max_target_seqs 10 -entrez_query "Sus scrofa [ORGN]"

#Convert the archive output to a tabular output
blast_formatter -archive UK1_Sus_scrofa_archive_070323.txt -outfmt 6 -out UK1_Sus_scrofa_tabular_070323.txt