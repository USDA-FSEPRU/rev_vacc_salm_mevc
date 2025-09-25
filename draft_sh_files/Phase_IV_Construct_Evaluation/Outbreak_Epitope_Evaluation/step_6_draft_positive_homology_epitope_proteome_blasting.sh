#!/bin/bash
#SBATCH --job-name=Step_6_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Pt1
#SBATCH -p short
#SBATCH --mem=4gb
#SBATCH --array=1-10
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 1:00:00
#SBATCH -o "Step_6_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Pt1.stdout.%j.%N"
#SBATCH -e "Step_6_Loc_Indpt_Epitopes_PulseNet_Proteomes_Blasting_Pt1.stderr.%j.%N"

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
proteome_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_bakta_protein_faas/*.faa | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
proteome_name="$(echo $proteome_file | cut -d "/" -f 8 | cut -d "." -f 1)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#proteome_file="$(ls /project/fsepru113/dbradshaw/reverse_vaccinology/Outbreak_Comparisons/PulseNet_bakta_protein_faas/*.faa | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#proteome_name="$(echo $proteome_file | cut -d "/" -f 8 | cut -d "." -f 1)"

source /project/fsepru113/dbradshaw/miniconda3/bin/activate blast

#Make a blast database
makeblastdb -in $proteome_file -title $proteome_name -dbtype prot -out blast_dbs/$proteome_name -parse_seqids

#Blast against database
blastp -query VaxiJen_Localization_Agnostic_Construct_Epitopes.fasta -db blast_dbs/$proteome_name -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle ssciname qcovhsp" -out epitope_blasting_results/"$proteome_name"_tabular_PulseNet.txt -evalue 10000 -max_target_seqs 5