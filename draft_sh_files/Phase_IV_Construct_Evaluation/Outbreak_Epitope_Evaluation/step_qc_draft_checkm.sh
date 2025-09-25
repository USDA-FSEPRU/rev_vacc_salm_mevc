#!/bin/bash
#SBATCH --job-name=Step_QC_Infantis_IB_PulseNet_19_23_Turkey_Human_USA_030124_checkM
#SBATCH -p short
#SBATCH --mem-per-cpu=2gb
#SBATCH -c 48
#SBATCH -o "Step_QC_Infantis_IB_PulseNet_19_23_Turkey_Human_USA_030124_checkM.stdout.%j.%N"
#SBATCH -e "Step_QC_Infantis_IB_PulseNet_19_23_Turkey_Human_USA_030124_checkM.stderr.%j.%N"

source /project/fsepru113/dbradshaw/miniconda3/bin/activate assembly_evaluation


checkm taxonomy_wf species "Salmonella enterica" /90daydata/fsepru113/dbradshaw/Senftenberg_Projects/Infantis_021524/all_assemblies /90daydata/fsepru113/dbradshaw/Senftenberg_Projects/Infantis_021524/checkm -t 48 --tab-table