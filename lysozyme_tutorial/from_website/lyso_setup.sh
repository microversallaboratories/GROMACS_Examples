#!/bin/bash

# Download the file from RCSB to the datapath directory
for id in ${GROMACS_PDB}
do
  wget https://files.rcsb.org/download/${id}.pdb
done

# Clean the file by removing water molecules
grep -v $PDB_REMOVE ${GROMACS_PDB}.pdb > ${GROMACS_PDB}_clean.pdb
echo "File cleaned of unwanted res by grep"

# Convert the file to a .gro file with pdb2gmx
gmx pdb2gmx -f ${GROMACS_PDB}_clean.pdb -o ${GROMACS_PDB}_processed.gro -water $GROMACS_WATERMODEL

# Edit the box by changing its dimensions
gmx editconf -f ${GROMACS_PDB}_processed.gro -o ${GROMACS_PDB}_newbox.gro -c -d 1.0 -bt cubic

# Solvate the box
gmx solvate -cp ${GROMACS_PDB}_newbox.gro -cs spc216.gro -o ${GROMACS_PDB}_solv.gro -p topol.top

# Rename temp.top*****
# mv temp.top* temp.top

# Download the requisite .mdp parameter file
wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp

# Generate the restraint file
gmx grompp -f ions.mdp -c ${GROMACS_PDB}_solv.gro -p topol.top -o ions.tpr

# Generate ions inside the box
gmx genion -s ions.tpr -o ${GROMACS_PDB}_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Download an input parameter file
wget http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp

# Assemble the binary input 
gmx grompp -f minim.mdp -c ${GROMACS_PDB}_solv_ions.gro -p topol.top -o em.tpr