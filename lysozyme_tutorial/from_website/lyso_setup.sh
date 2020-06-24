#!/bin/bash

# Download the file from RCSB
wget https://files.rcsb.org/view/1AKI.pdb

# Clean the file by removing water molecules
grep -v HOH 1aki.pdb > 1AKI_clean.pdb

# Convert the file to a .gro file with pdb2gmx
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce

# Edit the box by changing its dimensions
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic

# Solvate the box
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top

# Download the requisite .mdp parameter file
wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp

# Generate the restraint file
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr

# Generate ions inside the box
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Download an input parameter file
wget http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp

# Assemble the binary input 
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr