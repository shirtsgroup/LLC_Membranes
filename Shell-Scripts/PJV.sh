#!/bin/bash
# Script to build a structure to a specified number of layers using a specified
# monomer, then energy minimize the structure, and finally run a simulation
# based on inputs to the .mdp files

# Define Variables

ifv.sh

INTEGRATOR='steep'
NSTEPS_EM=50000
CUTOFF_EM='verlet'
NSTLIST=10

sed -i -e "s/INTEGRATOR/${INTEGRATOR}/g" em.mdp
sed -i -e "s/NSTEPS_EM/${NSTEPS_EM}/g" em.mdp
sed -i -e "s/CUTOFF_EM/${CUTOFF_EM}/g" em.mdp
sed -i -e "s/NSTLIST_EM/${NSTLIST}/g" em.mdp

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

INPUT_FILE=$DIR/../Structure-Files/monomer4.pdb
LAYERS=20    # number of layers wanted in the structure
Z_BOX_VECTOR=$((LAYERS+5)) # kind of arbitrary but should work
XVECT=8.0
YVECT=8.0
INCREMENT=0.1

# Build Structure Based on user-defined inputs

python $DIR/../Structure_Builder/Structure_Builder_for_Bash.py -i $INPUT_FILE -l $LAYERS >> initial.gro

# Put the structure in a box

gmx editconf -f initial.gro -o box.gro -c -bt triclinic -box $XVECT $YVECT $Z_BOX_VECTOR -angles 90 90 120

# Prepare input file (.tpr) for energy minimization

gmx grompp -f em.mdp -c box.gro -p NaPore.top -o box_em.tpr

# Run Energy Minimization

gmx mdrun -v -deffnm box_em

# Extract Potential Energy from log file
ENERGY=$(cat box_em.log | grep 'Potential Energy' | awk '{print substr($0,21,5)}')

# If the potential energy is positive, then the box vector is incremented by a fixed amount until the energy comes out negative

while [ $(echo " $ENERGY > 0" | bc) -eq 1 ]; do
        XVECT=$(echo "$XVECT + $INCREMENT" | bc -l)
        YVECT=$(echo "$YVECT + $INCREMENT" | bc -l)
        gmx editconf -f initial.gro -o box.gro -c -bt triclinic -box $XVECT $YVECT $Z_BOX_VECTOR -angles 90 90 120
        gmx grompp -f em.mdp -c box.gro -p NaPore.top -o box_em.tpr
        gmx mdrun -v -deffnm box_em
        ENERGY=$(cat box_em.log | grep 'Potential Energy' | awk '{print substr($0,21,5)}')
done

# prepare input file for simulation (just letting it wiggle around)

gmx grompp -f wiggle.mdp -c box_em.gro -p NaPore.top -o wiggle.tpr

# remove unecessary files
find . -type f -name 'box_em'\* -exec rm {} \;
find . -type f -name '#'\* -exec rm {} \;
#rm initial.gro  
rm box.gro 
