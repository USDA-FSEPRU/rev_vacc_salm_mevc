#!/bin/bash
#SBATCH --job-name=Step_4_PulseNet_Enteritidis_skesa_assembly
#SBATCH -p short
#SBATCH --mem-per-cpu=12gb
#SBATCH -c 4
#SBATCH --array=1-3
#SBATCH -o "Step_4_PulseNet_Enteritidis_skesa_assembly.stdout.%j.%N"
#SBATCH -e "Step_4_PulseNet_Enteritidis_skesa_assembly.stderr.%j.%N"

read_id="$(ls reads/fastp/ | cut -d "." -f 1 | sort | uniq | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"


source /project/fsepru113/dbradshaw/miniconda3/bin/activate skesa

mkdir skesa_assembly/

skesa --reads reads/fastp/"$read_id".fastp_1.fastq,reads/fastp/"$read_id".fastp_2.fastq --cores 4 --memory 48 > skesa_assembly/"$read_id".skesa.fna