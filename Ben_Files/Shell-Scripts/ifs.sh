#!/bin/bash
# Script to copy over all necessary files for simulation of solvated system
# Edit files copied in as necessary

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

for file in wiggle.mdp wiggle_solv.mdp NaPore.top NaPore_water.top vdwradii.dat em.mdp; do
	cp $DIR/Files/$file $PWD

done
