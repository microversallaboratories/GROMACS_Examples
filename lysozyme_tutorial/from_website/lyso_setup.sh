#!/bin/bash

#This script needs to be edited for each run.

# Define the timestamp for the new simulation, to identify the simulation folder and files
JLN_SIM_TIMESTAMP=$(date +"%Y_%m_%d_%H_%M_%S")

# Define the simulation type (prefix)
JLN_SIM_TYPE='mdsim'

# Define PDB Filename, paths, & GROMACS Pameters
JLN_SIM_PDBNAME='1aki'			# Should eventually take this in at runtime!
JLN_DATAPATH=/Box/extracurriculars/research/SURE_S2020_fileshare/sure_data		# Datafile path
JLN_SIM_PATH=$JLN_DATAPATH/$JLN_SIM_PDBNAME/"${JLN_SIM_TYPE}_${JLN_SIM_TIMESTAMP}"		# Simulation-specific datafile path
JLN_ANAPATH=$JLN_SIM_PATH/analysis		# Analysis filepath

# Make a directory for the current simulation (ANAPATH creates all three target directories)
mkdir -p /Users/jacobnorth/"$JLN_ANAPATH"

# Gromacs parameters
GROMACS_PDB=$1				# Take in first cmdlnarg
PDB_REMOVE="HOH"
GROMACS_FORCEFIELD="gromos53a6"		# Define ff name
GROMACS_WATERMODEL="spc"			# Define water model name
#GROMACS_BOXTYPE="dodecahedron"		# Def boxtype
GROMACS_BOXTYPE="cubic"
GROMACS_BOXORIENTATION="1.5"		# Def box orientation
GROMACS_BOXSIZE="5.0"				# Def boxsize
GROMACS_BOXCENTER="2.5"				# Def boxcenter

#Setup GROMACS Job. Probably not necessary to edit past this point.
if [ -z "$JLN_SIM_PATH/$GROMACS_PDB" ]; then
	echo "USAGE: ./setup_GROMACS_job.sh pdb_filename"
	echo "Do NOT include the .pdb extension in the file name."
	exit
fi

# Download the file from RCSB to the datapath directory
wget -O $JLN_SIM_PATH/$GROMACS_PDB.pdb https://files.rcsb.org/view/$GROMACS_PDB.pdb

# Clean the file by removing water molecules
grep -v $PDB_REMOVE $JLN_SIM_PATH/${GROMACS_PDB}.pdb > $JLN_SIM_PATH/${GROMACS_PDB}_clean.pdb

# Convert the file to a .gro file with pdb2gmx
gmx pdb2gmx -f $JLN_SIM_PATH/${GROMACS_PDB}_clean.pdb -o $JLN_SIM_PATH/${GROMACS_PDB}_processed.gro -water $GROMACS_WATERMODEL

# Edit the box by changing its dimensions
gmx editconf -f $JLN_SIM_PATH/${GROMACS_PDB}_processed.gro -o $JLN_SIM_PATH/${GROMACS_PDB}_newbox.gro -c -d 1.0 -bt cubic

# Solvate the box
gmx solvate -cp $JLN_SIM_PATH/${GROMACS_PDB}_newbox.gro -cs $JLN_SIM_PATH/spc216.gro -o $JLN_SIM_PATH/${GROMACS_PDB}_solv.gro -p $JLN_SIM_PATH/topol.top

# Download the requisite .mdp parameter file
wget -O $JLN_SIM_PATH/ions.mdp http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp

# Generate the restraint file
gmx grompp -f $JLN_SIM_PATH/ions.mdp -c $JLN_SIM_PATH/${GROMACS_PDB}_solv.gro -p $JLN_SIM_PATH/${GROMACS_PDB}_topol.top -o $JLN_SIM_PATH/ions.tpr

# Generate ions inside the box
gmx genion -s $JLN_SIM_PATH/ions.tpr -o $JLN_SIM_PATH/${GROMACS_PDB}_solv_ions.gro -p $JLN_SIM_PATH/${GROMACS_PDB}_topol.top -pname NA -nname CL -neutral

# Download an input parameter file
wget -O $JLN_SIM_PATH/minim.mdp http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp

# Assemble the binary input 
gmx grompp -f $JLN_SIM_PATH/minim.mdp -c $JLN_SIM_PATH/${GROMACS_PDB}_solv_ions.gro -p $JLN_SIM_PATH/topol.top -o $JLN_SIM_PATH/em.tpr