#!/bin/bash
# Script to build a structure to a specified number of layers using a specified
# monomer, then energy minimize the structure, run a simulation, then readjust the box size to give a desired number of nm of water between layers, energy minimize again, and run a simulation.

# based on inputs to the .mdp files

# Define Variables

ifs.sh

INPUT_FILE=/home/bcoscia/PycharmProjects/LLC_Structure_Builder/Monomer_Configurations/monomer4.pdb
LAYERS=20    # number of layers wanted in the structure
Z_BOX_VECTOR=$((LAYERS+5)) # kind of arbitrary but should work
WATER_LAYER=6  # nm of water between layers in the z direction
# Starting x and y box vectors:
XVECT1=8.5  # for first energy minimization
YVECT1=8.5
XVECT2=8.5  # for second energy minimization
YVECT2=8.5
INCREMENT=0.1  # How much each box vector will be increased in the event of a positive energy found through energy minimzation

# Build Structure Based on user-defined inputs

python /home/bcoscia/PycharmProjects/LLC_Structure_Builder/Structure_Builder_for_Bash.py -i $INPUT_FILE -l $LAYERS >> initial.gro

# Put the structure in a box

gmx editconf -f initial.gro -o box.gro -c -bt triclinic -box $XVECT1 $YVECT1 $Z_BOX_VECTOR -angles 90 90 120

# Prepare input file (.tpr) for energy minimization

gmx grompp -f em.mdp -c box.gro -p NaPore.top -o box_em.tpr

# Run Energy Minimization

gmx mdrun -v -deffnm box_em

# Extract Potential Energy from log file
ENERGY1=$(cat box_em.log | grep 'Potential Energy' | awk '{print substr($0,21,5}')

# If the potential energy is positive, then the box vector is incremented by a fixed amount until the energy comes out negative

while [ $(echo " $ENERGY1 > 0" | bc) -eq 1 ]; do
        XVECT1=$(echo "$XVECT1 + $INCREMENT" | bc -l)
        YVECT1=$(echo "$YVECT1 + $INCREMENT" | bc -l)
        gmx editconf -f initial.gro -o box.gro -c -bt triclinic -box $XVECT1 $YVECT1 $THICKNESS -angles 90 90 120
        gmx grompp -f em.mdp -c box.gro -p NaPore.top -o box_em.tpr
        gmx mdrun -v -deffnm box_em
	ENERGY1=$(cat box_em.log | grep 'Potential Energy' | awk '{print substr($0,21,5)}')
done

# prepare input file for simulation (just letting it wiggle around)

gmx grompp -f wiggle.mdp -c box_em.gro -p NaPore.top -o wiggle_vac.tpr

# run simulation for amount of time specified in wiggle.mdp

gmx mdrun  -v -deffnm wiggle_vac

INPUT_VAC_FILE=wiggle_vac.gro

THICKNESS=$(python /home/bcoscia/PycharmProjects/Membrane_Dimensions/Thickness_Bash.py -i $INPUT_VAC_FILE -w $WATER_LAYER)

# Edit box to correct dimensions (i.e. so there is room for the specified amount of water)

gmx editconf -f wiggle_vac.gro -o new_box.gro -c -bt triclinic -box 8.5 8.5 $THICKNESS -angles 90 90 60

# energy minimize in this new box since everything is shifted around

gmx grompp -f em.mdp -c new_box.gro -p NaPore.top -o em_new_box.tpr

gmx mdrun -v -deffnm em_new_box

# Extract Potential energy from that minimization run. A positive value indicates that the box vectors are too small

ENERGY2=$(cat em_new_box.log | grep 'Potential Energy' | awk '{print substr($0,21,5)}')
# if the potential energy is positive, then the box vector is increased by a defined increment until the energy is negative after energy minimzation

while [ $(echo "  $ENERGY2 > 0" | bc) -eq 1 ]; do
        XVECT2=$(echo "$XVECT2 + $INCREMENT" | bc -l)
        YVECT2=$(echo "$YVECT2 + $INCREMENT" | bc -l)
        gmx editconf -f wiggle_vac.gro -o new_box.gro -c -bt triclinic -box $XVECT2 $YVECT2 $THICKNESS -angles 90 90 60
        gmx grompp -f em.mdp -c new_box.gro -p NaPore.top -o em_new_box.tpr
        gmx mdrun -v -deffnm em_new_box
	ENERGY2=$(cat em_new_box.log | grep 'Potential Energy' | awk '{print substr($0,21,5)}')
done

# run simulation in new energy minimized box (Correct_box)

gmx grompp -f wiggle.mdp -c em_new_box.gro -p NaPore.top -o Correct_box.tpr

gmx mdrun -v -deffnm Correct_box

# Now solvate the box

gmx solvate -cp Correct_box.gro -cs spc216.gro -p NaPore_water.top -o water.gro

# Now run a simulation with water in the box
# first energy minimize

gmx grompp -f em.mdp -c water.gro -p NaPore_water.top -o water_em.tpr

gmx mdrun -v -deffnm water_em

# Now run the actual simulation

gmx grompp -f wiggle_solv.mdp -c water_em.gro -p NaPore_water.top -o water_wiggle.tpr

gmx mdrun -v -deffnm water_wiggle

gmx trjconv -f water_wiggle.trr -s water_wiggle.tpr -o traj.gro

# remove unecessary files
rm initial.gro water.gro new_box.gro #remove now unnecessary file
rm box.gro # remove unecessary file
find . -type f -name 'Correct_box'\* -exec rm {} \;
find . -type f -name 'water_em'\* -exec rm {} \;
find . -type f -name 'em_new_box'\* -exec rm {} \;
find . -type f -name 'box_em'\* -exec rm {} \;
find . -type f -name '#'\* -exec rm {} \;
find . -type f -name 'wiggle_vac'\* -exec rm {} \;
# Directory for output files that aren't needed that often
mkdir Output
mv mdout.mdp water_wiggle.cpt water_wiggle.edr water_wiggle.log water_wiggle.tpr Output
