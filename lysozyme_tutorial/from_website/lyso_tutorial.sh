#!/bin/bash

#This script needs to be edited for each run.

#%%%%%%%%%%%%%%%%%%%%<CONFIGURATION>%%%%%%%%%%%%%%%%%%%#

# Define the timestamp for the new simulation, to identify the simulation folder and files
JLN_SIM_TIMESTAMP=$(date +"%Y_%m_%d_%H_%M_%S")

# Define the simulation type (prefix)
JLN_SIM_TYPE='mdsim'

# Define PDB Filename, paths, & GROMACS Pameters
JLN_SIM_PDBNAME='1aki'			# Should eventually take this in at runtime!
JLN_DATAPATH=/Users/jacobnorth/Box/extracurriculars/research/SURE_S2020_fileshare/sure_data		# Datafile path
JLN_SIM_PATH=$JLN_DATAPATH/$JLN_SIM_PDBNAME/"${JLN_SIM_TYPE}_${JLN_SIM_TIMESTAMP}"		# Simulation-specific datafile path
JLN_ANAPATH=$JLN_SIM_PATH/analysis		# Analysis filepath

# Make a directory for the current simulation (ANAPATH creates all three target directories)
mkdir -p /Users/jacobnorth/"$JLN_ANAPATH"

echo "Simulation filepath is:"
echo ${JLN_SIM_PATH}

#%%%%%%%%%%%%%%%%%%%</CONFIGURATION>%%%%%%%%%%%%%%%%%%%#

#%%%%%%%%%%%%%%%%%%%<GROMACS CONF>%%%%%%%%%%%%%%%%%%%%%#

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

# Setup GROMACS Job. Probably not necessary to edit past this point.
if [ -z "$JLN_SIM_PATH/$GROMACS_PDB" ]; then
	echo "USAGE: ./setup_GROMACS_job.sh pdb_filename"
	echo "Do NOT include the .pdb extension in the file name."
	exit
fi

#%%%%%%%%%%%%%%%%%%%</GROMACS CONF>%%%%%%%%%%%%%%%%%%%%#

# Then, just call all those scripts in order :)

# CLEANUP 
bash lyso_cleanup.sh			# Clean up previous residual files

# SIMULATION SETUP
bash lyso_setup.sh				# Setup the simulation
bash lyso_emin.sh				# Minimize the energy of the system
bash lyso_nvtequilib.sh			# Equilibrate wrt NVT
bash lyso_nptequilib.sh			# Equilibrate wrt NPT

# PRODUCTION MD
bash lyso_prodmd.sh				# Run production MD
bash lyso_postpmd_analysis.sh	# Analyze post production MD

# RELOCATE FILES TO THE TIMESTAMPED DIRECTORY
bash lyso_relocate_sim_files.sh	# Relocate files