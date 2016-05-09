#!/bin/bash
# Script to copy over all necessary files for simulation of system in a vacuum
# Edit files copied in as necessary

for file in wiggle.mdp NaPore.top em.mdp; do
	cp /home/bcoscia/Documents/Gromacs/Sim_Scripts/Files/$file $PWD

done
