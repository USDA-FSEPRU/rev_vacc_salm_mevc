#!/bin/bash
#SBATCH --job-name=Typhimurium_proteome_blasting
#SBATCH -p short
#SBATCH --mem=4gb
#SBATCH --array=1-2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 1:00:00
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"


proteome_file="$(ls protein_data/*.faa | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
proteome_name="$(echo $proteome_file | cut -d "/" -f 2 | cut -d "_" -f 1,2)"


source /project/fsepru113/dbradshaw/miniconda3/bin/activate blast

#Make a blast database
makeblastdb -in $proteome_file -title $proteome_name -dbtype prot -out blast_dbs/$proteome_name -parse_seqids

#Blast against database
blastp -query Pre_Positive_Homology_UK1_Epitopes.fasta -db blast_dbs/$proteome_name -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle ssciname" -out epitope_blasting_results/"$proteome_name"_tabular_070723.txt -evalue 10000 -max_target_seqs 5