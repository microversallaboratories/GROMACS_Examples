#!/bin/bash


# Invoke the energy minimization
gmx mdrun -v -deffnm em

# Plot energy chart for eminim 
gmx energy -f $GROMACS_JLN_CURR_DATAPATH/1aki/mdsim_2020_6_24_15_44_SHELLSCRIPT/em.edr -o potential.xvg