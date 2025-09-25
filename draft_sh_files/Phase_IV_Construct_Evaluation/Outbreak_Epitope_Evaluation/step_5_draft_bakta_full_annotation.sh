#!/bin/bash
#SBATCH --job-name=Step_5_PulseNet_All_Bakta_Annotation
#SBATCH -p short
#SBATCH --mem-per-cpu=6gb
#SBATCH --array=1-10
#SBATCH -c 4
#SBATCH -o "Step_5_PulseNet_All_Bakta_Annotation.stdout.%j.%N"
#SBATCH -e "Step_5_PulseNet_All_Bakta_Annotation.stderr.%j.%N"

#Activate the environment
source /project/fsepru113/dbradshaw/miniconda3/bin/activate bakta

#Create a variable to sequentially store assembly file and sample id names (first 10k set)
assembly_file="$(ls /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/final_assemblies/ | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"
sample_id="$(echo $assembly_file| cut -d "." -f 1)"

#Create a variable to sequentially store assembly file and sample id names (next 10k sets)
#assembly_file="$(ls /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/final_assemblies/ | head -n $((SLURM_ARRAY_TASK_ID+10000)) | tail -n 1)"
#sample_id="$(echo $assembly_file| cut -d "." -f 1)"


#Example of database download
#bakta_db download --output /project/fsepru113/dbradshaw/databases/bakta_db_light --type light

#Run bakta to annotate each subseted assembly of pESI-like contigs using a complete pESI reference
bakta --db /project/fsepru113/dbradshaw/databases/bakta_db_full/ -o bakta_full_annotation/$sample_id.bakta --keep-contig-headers -t 4 --genus Salmonella --species enterica --verbose --force --prefix $sample_id /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/final_assemblies/$assembly_file