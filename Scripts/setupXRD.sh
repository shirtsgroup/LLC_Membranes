#!/bin/bash

set -e

IMAGES=1
LAYERS=20
BUILD_MON='NAcarb11Vd'
Z=15

while getopts "I:l:b:z:" opt; do
    case $opt in
    I) IMAGES=$OPTARG;;
    l) LAYERS=$OPTARG;;
    b) BUILD_MON=$OPTARG;;
    z) Z=$OPTARG;;
    esac
done

# Set up system and simulation files
build.py -l ${LAYERS}
gmx editconf -f initial.gro -o box.gro -c -bt triclinic -box 9 9 ${Z}
input.py -b ${BUILD_MON} -e 'npt' -p 'isotropic' -l 50 -r 'on' -c 'box.gro'

# equilibrate
equil.sh

# duplicate periodically
Periodic_Images.py -f 0.gro -o '${IMAGES}images.gro' -i ${IMAGES}

gmx grompp -f em.mdp -p topol2.top -c periodic.gro -o periodic
gmx mdrun -v -deffnm periodic

gmx grompp -f npt.mdp -p topol2.top -c periodic.gro -o wiggle
gmx mdrun -v -deffnm wiggle

