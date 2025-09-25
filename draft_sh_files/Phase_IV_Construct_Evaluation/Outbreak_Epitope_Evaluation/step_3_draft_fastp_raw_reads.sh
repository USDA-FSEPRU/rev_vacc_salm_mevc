#!/bin/bash
#SBATCH --job-name=Step_3_PulseNet_Enteritidis_fastp_raw_reads
#SBATCH -p short
#SBATCH --mem-per-cpu=4gb
#SBATCH -c 4
#SBATCH --array=1-2
#SBATCH -o "Step_3_PulseNet_Enteritidis_fastp_raw_reads.stdout.%j.%N"
#SBATCH -e "Step_3_PulseNet_Enteritidis_fastp_raw_reads.stderr.%j.%N"

read_id="$(ls reads/raw/ | cut -d "_" -f 1 | sort | uniq | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)"

source /project/fsepru113/dbradshaw/miniconda3/bin/activate polish

fastp --in1 reads/raw/"$read_id"_1.fastq --in2 reads/raw/"$read_id"_2.fastq --out1 reads/fastp/"$read_id".fastp_1.fastq --out2 reads/fastp/"$read_id".fastp_2.fastq --unpaired1 reads/unpaired/"$read_id".fastp.unpaired_1.fastq --unpaired2 reads/unpaired/"$read_id"_fastp.unpaired_2.fastq