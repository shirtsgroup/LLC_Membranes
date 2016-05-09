#!/bin/bash
# Script to copy over all necessary files for simulation of solvated system
# Edit files copied in as necessary

for file in wiggle.mdp wiggle_solv.mdp NaPore.top NaPore_water.top vdwradii.dat em.mdp; do
	cp /home/bcoscia/Documents/Gromacs/Sim_Scripts/Files/$file $PWD

done
