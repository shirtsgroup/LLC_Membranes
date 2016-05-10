#!/bin/bash
# Script to build a structure to a specified number of layers using a specified
# monomer, then energy minimize the structure, and finally run a simulation
# based on inputs to the .mdp files

# Copy in necessary files
ifv.sh

# Define Variables

# Choose which monomer to build with
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"  # Directory where this script is located
MONOMER='monomer4.pdb'  # Structure file to be used
INPUT_FILE=$DIR/../Structure-Files/$MONOMER  # Full path to file where monomer structure is located

# Energy minimization parameters:
INTEGRATOR='steep'  # Integrator for energy minimization
NSTEPS_EM=50000  # Maximum number of steps to take for energy minimization
CUTOFF_EM='verlet'  # Cut-off Scheme
NSTLIST=10  # Neighborlist - changed automatically by gromacs unless it is set equal to 1

# Membrane Dimensions
NO_MONOMERS=6  # Number of monomers in 1 layer
RADIUS=3  # Initial pore radius, angstroms
PORE2PORE=40  # Pore-to-Pore distance, angstroms
NOPORES=4  # Number of pores to be built
DBWL=10  # Distance Between Layers
LAYERS=20    # Number of layers wanted in the structure

# Box Vector Parameters
Z_BOX_VECTOR=$((LAYERS+5)) # kind of arbitrary but should work
XVECT=8.0  # Box vector in the x direction
YVECT=8.0  # Box vector in the y direction (not this will be multiplied by sin(120) for a monoclinic cell of 120 degrees
INCREMENT=0.1  # Increment to increase the box vector by if there is a LINCS error

sed -i -e "s/INTEGRATOR/${INTEGRATOR}/g" em.mdp
sed -i -e "s/NSTEPS_EM/${NSTEPS_EM}/g" em.mdp
sed -i -e "s/CUTOFF_EM/${CUTOFF_EM}/g" em.mdp
sed -i -e "s/NSTLIST_EM/${NSTLIST}/g" em.mdp



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
