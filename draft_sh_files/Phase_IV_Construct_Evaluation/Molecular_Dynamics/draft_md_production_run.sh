#!/bin/bash
#SBATCH --job-name=Loc_Indpdt_TLR5_Complex_md_prod_run
#SBATCH -p long
#SBATCH --mem=144gb
#SBATCH -N 1
#SBATCH -n 72
#SBATCH -o "Loc_Indpdt_TLR5_Complex_md_prod_run.stdout.%j.%N"
#SBATCH -e "Loc_Indpdt_TLR5_Complex_md_prod_run.stderr.%j.%N"

#Run with sbatch Loc_Indpdt_TLR5_Complex_md_prod_run.sh --no-requeue

#Activate gromacs
source /project/fsepru113/dbradshaw/gromacs-2023.3/bin/GMXRC

#Run a full simulation - First Time
gmx mdrun -v -deffnm md

#Run a continuation of a simulation using a checkpoint
#gmx mdrun -v -deffnm md -cpi md.cpt