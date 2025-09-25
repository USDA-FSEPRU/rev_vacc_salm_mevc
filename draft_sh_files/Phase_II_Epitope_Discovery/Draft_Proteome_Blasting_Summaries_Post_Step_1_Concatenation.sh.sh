#!/bin/bash
#SBATCH --job-name=PulseNet_combining_summaries
#SBATCH -p short
#SBATCH --mem-per-cpu=12gb
#SBATCH -c 1
#SBATCH -o "PulseNet_combining_summaries.stdout.%j.%N"
#SBATCH -e "PulseNet_combining_summaries.stderr.%j.%N"

head -n 1 SRR10000867_proteome_epitope_homology.txt > PulseNet_proteomes_epitope_homology.txt; tail -n +2 -q *.txt >> PulseNet_proteomes_epitope_homology.txt
