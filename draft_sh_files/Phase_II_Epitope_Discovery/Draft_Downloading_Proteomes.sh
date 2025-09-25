#!/bin/bash
#SBATCH --job-name=Enteritids_proteomes
#SBATCH -p medium
#SBATCH --mem=12gb
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"

source /project/fsepru113/dbradshaw/miniconda3/bin/activate edirect

esearch -db assembly -query txid149539[organism:exp] | esummary   | xtract -pattern DocumentSummary -element FtpPath_GenBank   | while read -r url ; do      path=$(echo $url | perl -pe 's/(GC[FA]_\d+.*)/\1\/\1_protein.faa.gz/g') ;  wget -q --show-progress "$path" -P protein_data ;    done