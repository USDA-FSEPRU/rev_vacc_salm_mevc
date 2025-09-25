#!/bin/bash
#SBATCH --job-name=Loc_Indpdt_TLR5_Complex_md_prod_run
#SBATCH -p long
#SBATCH --mem=144gb
#SBATCH -N 1
#SBATCH -n 72
#SBATCH -o "Loc_Indpdt_TLR5_Complex_md_prod_run.stdout.%j.%N"
#SBATCH -e "Loc_Indpdt_TLR5_Complex_md_prod_run.stderr.%j.%N"

#Activate gromacs
source /project/fsepru113/dbradshaw/gromacs-2023.3/bin/GMXRC

##Generate graphics

#Generate a graphic for RMSD - < 5 min (10 ns sim)
#https://www.compchems.com/what-is-the-rmsd-and-how-to-compute-it-with-gromacs/
#-s is the reference geometry
#-f is the trajectory file
#-tu ns will change the unites to ns instead of ps
printf "4 4" | gmx rms -f final.xtc -s md_50.tpr -o md_50_rmsd.xvg
#Choose 4 to select backbone group for least squares fit and 4 to choose backbone group for RMSD calculation

#Generate a graphic for radius of gyration - < 3 min (10 ns sim)
#https://tutorials.gromacs.org/docs/md-intro-tutorial.html
printf "1" | gmx gyrate -f final.xtc -s md_50.tpr -o md_50_gyrate.xvg
#Choose 1  to select protein group

#Generate a graphic for RMSF (Root Mean Square Fluctuation) - < 3 min (10 ns sim)
#https://www.compchems.com/how-to-compute-the-rmsf-using-gromacs/#what-is-the-rmsf
#From https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html - An area of the structure with high RMSF values frequently diverges from the average, indicating high mobility. When RMSF analysis is carried out on proteins, it is typically restricted to backbone or alpha-carbon atoms; these are more characteristic of conformational changes than the more flexible side-chains.
#From https://gromacs.bioexcel.eu/t/rmsf-analysis-of-protein-ligand-md-simulation/4919/5 - I would also recommend performing RMSF analysis only on the Cα or backbone groups, not the whole protein. You will get a lot of irrelevant motions (like sidechains) if you consider everything.
# -f is the trajectory file
# -s is the reference structure
# -res is to compute the RMSF for each Residue
# -b and -e allow you to choose a specific time frame in ps
printf "3" | gmx rmsf -f final.xtc -s md_50.tpr -o md_50_rmsf.xvg -res
#Choose 3 to select the C-alpha group