#!/bin/bash
#SBATCH --job-name=Step_2_PulseNet_Enteritidis_SRAs_to_Reads
#SBATCH -p short
#SBATCH --mem-per-cpu=1gb
#SBATCH -c 6
#SBATCH --array=1-2
#SBATCH -o "Step_2_PulseNet_Enteritidis_SRAs_to_Reads.stdout.%j.%N"
#SBATCH -e "Step_2_PulseNet_Enteritidis_SRAs_to_Reads.stderr.%j.%N"

sra="$(ls sras/*.sra* | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"

#Load SRA Toolkit
module load sratoolkit

#Split all sra files to paired reads
fasterq-dump -S $sra -O reads/raw --verbose -e 6