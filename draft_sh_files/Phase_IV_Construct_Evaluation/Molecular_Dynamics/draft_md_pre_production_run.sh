#!/bin/bash
#SBATCH --job-name=Loc_Dpdt_TLR2_md_prod_run
#SBATCH -p short
#SBATCH --mem=144gb
#SBATCH -N 1
#SBATCH -n 72
#SBATCH -o "Loc_Dpdt_TLR2_md_prod_run.stdout.%j.%N"
#SBATCH -e "Loc_Dpdt_TLR2_md_prod_run.stderr.%j.%N"

#Activate gromacs
source /project/fsepru113/dbradshaw/gromacs-2023.3/bin/GMXRC

#Delete water and ligands from pdb file - < 1 min MM; < 1 min DM
grep -v HETATM TLR5_Homodimer_Any_Loc_Refined_Balanced_ClusPro_Model.000.06.pdb > TLR5_Any_clean.pdb

#Generate a topology file - < 1 min MM; < 1 min DM
printf "15" | gmx pdb2gmx -f TLR5_Any_clean.pdb -o TLR5_Any_clean.gro -water spce -ignh
#Choose option 15 for OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)
#Generates a posre file for each chain, a .gro file, a *.itp file for each chain, and a topol.top file

#Define a box - < 1 min MM; < 1 min DM
gmx editconf -f TLR5_Any_clean.gro -o TLR5_Any_box.gro -c -d 1.0 -bt cubic
#Generates the TLR5_Any_box.gro file

#Solvate the system - < 1 min MM; < 1 min DM
gmx solvate -cp TLR5_Any_box.gro -cs spc216.gro -o TLR5_Any_solv.gro -p topol.top
#Generates the TLR5_Any_solv.gro file

##Neutralize the system

#Assemble the ions.tpr file  - < 1 min MM; < 1 min DM
gmx grompp -f ions.mdp -c TLR5_Any_solv.gro -p topol.top -o ions.tpr
#Generates the ions.tpr file and another topol.top file

#Genion module  - < 1 min MM; < 1 min DM
printf "13" | gmx genion -s ions.tpr -o TLR5_Any_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
#Choose option 13 - SOL to ensure that ions only appear in the solvent
#Generates the TLR5_Any_ions.gro file and a new topol.top file

##Energy minimization

#Generate tpr file - < 1 min MM; < 1 min DM
gmx grompp -f minim.mdp -c TLR5_Any_ions.gro -p topol.top -o em.tpr
#Generates the em.tpr file

#Run simulation - < 3 hrs MM; < 1.5 hrs DM
gmx mdrun -v -deffnm em
#Generates em.log, em.gro, em.edr, and em.trr files

#Generate a plot of potential over time - < 1 min MM; < 1 min DM
printf "10 0" | gmx energy -f em.edr -o potential.xvg
#Choose 10 for Potential Energy, then 0 to end the selection
#Generates the potential.xvg file

##NVT equilibration

#Generate a tpr file - < 1 min MM; < 1 min DM
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
#Generates teh nvt.tpr file

#Run simulation - <  MM; < 1 hr DM
gmx mdrun -v -deffnm nvt
#Generates the nvt.edr, nvt.gro, nvt.log, nvt.trr, nvt_prev.cpt, and nvt.cpt files

#Generate a graphic for GROMACS energy - < 1 min MM; < 1 min DM
printf "16 0" | gmx energy -f nvt.edr -o nvt_temperature.xvg
#Choose 16 for Temperature, then 0 to end the selection
#Generates the nvt_temperature.xvg file

#Generate a graphic for RMSD - < 1 min MM; < 1 min DM
printf "4 4" | gmx rms -f nvt.trr -s nvt.tpr -o nvt_rmsd.xvg
#Choose 4 to select backbone group for least squares fit and 4 to choose backbone group for RMSD calculation
#Generates the nvt_rmsd.xvg file

##NPT equilibration

#Generate a tpr file - < 1 min MM; < 1 min DM
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
#Generates the npt.tpr file

#Run simulation 
gmx mdrun -v -deffnm npt
#Generates the npt.edr, npt.gro, npt.log, npt.trr, npt_prev.cpt, and npt.cpt files

#Generate a graphic for GROMACS energy - < 1 min
printf "18 0" | gmx energy -f npt.edr -o npt_pressure.xvg
#Choose 18 for Pressure, then 0 to end the selection
#Generates the npt_pressure.xvg file

#Generate a graphic for RMSD - < 1 min
printf "4 4" | gmx rms -f npt.trr -s npt.tpr -o npt_rmsd.xvg
#Choose 4 to select backbone group for least squares fit and 4 to choose backbone group for RMSD calculation
#Generates the npt_rmsd.xvg file


##Production Run

#Generate a tpr file - < 1 min
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
#Generates the md.tpr file