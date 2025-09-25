#!/bin/bash
#SBATCH --job-name=Step_1_PulseNet_Enteritidis_SRRs
#SBATCH -p short
#SBATCH --mem-per-cpu=2gb
#SBATCH -c 20
#SBATCH -o "Step_1_PulseNet_Enteritidis_SRRs.stdout.%j.%N"
#SBATCH -e "Step_1_PulseNet_Enteritidis_SRRs.stderr.%j.%N"

#Load SRA Toolkit
module load sratoolkit

#Get all the sras in your list
prefetch --option-file PulseNet_Enteritidis_SRRs.txt --verbose

#Everything is downloaded into its own folder so make a new directory
mkdir sras

#Then move all sra files to the new folder 
cat PulseNet_Enteritidis_SRRs.txt | while read srr; do mv /90daydata/fsepru113/dbradshaw/reverse_vacc_salm_mevc/Outbreak_Comparisons/"$srr"/* sras/; done