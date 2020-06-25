#!/bin/bash

# Define MD datapath and MD anapath
DATAPATH=/Users/jacobnorth/Box/extracurriculars/research/SURE_S2020_fileshare/sure_data
SIMPATH=$DATAPATH/
ANAPATH=/Users/jacobnorth/Box/extracurriculars/research/SURE_S2020_fileshare/sure_data

# Run an nvt equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr