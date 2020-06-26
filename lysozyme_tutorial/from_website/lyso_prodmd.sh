#!/bin/bash

# Download the .mdp for the NVT equilib
wget http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp

# Pre-process the input files, prediec PME load
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

# Run the NVT equilibration:
gmx mdrun -deffnm md_0_1	# ON CPU
# gmx mdrun -deffnm md_0_1 -nb gpu	# ON ONE GPU
# ON MULTIPLE GPUS???

# Create plot of pressure, ensure T is stable
gmx energy -f npt.edr -o pressure.xvg

# Create plot of density
gmx energy -f npt.edr -o density.xvg