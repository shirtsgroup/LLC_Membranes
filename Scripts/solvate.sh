#!/usr/bin/env bash

WATER_LAYER_THICKNESS=3  # nm of water between membrane layers
GRO="wiggle.gro"  # name of .gro file to solvate
SOLVENT_CONFIGURATION="spc216.gro"
BUILD_MONOMER="NAcarb11V"
TOP="topol.top"
OUTPUT="wiggle"

while getopts "t:g:s:b:p:o:" opt; do
    case $opt in
    t) WATER_LAYER_THICKNESS=$OPTARG;;
    g) GRO=$OPTARG;;
    s) SOLVENT_CONFIGURATION=$OPTARG;;
    b) BUILD_MONOMER=$OPTARG;;
    p) TOP=$OPTARG;;
    o) OUTPUT=$OPTARG;;
    esac
done

input.py -b ${BUILD_MONOMER} -c ${GRO} -l 50 --solvate

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # the directory where this script is located.
cp ${DIR}/../top/vdwradii.dat .

THICKNESS=$(python -c "from llclib import physical; print physical.thickness('${GRO}')[0]")
ZBOX=`echo "${THICKNESS} + ${WATER_LAYER_THICKNESS}" | bc -l`

XBOX=$(tail -n 1 ${GRO} | head -c 10)  # get the last line of the gro file to get the box dimensions
# for simplicity, I'm going to assume XBOX only has one component and YBOX=XBOX.

# Now Gromacs can do the rest
gmx editconf -f ${GRO} -o extended.gro -bt triclinic -c -box ${XBOX} ${XBOX} ${ZBOX} -angles 90 90 120
gmx solvate -cp extended.gro -o solv.gro -cs ${SOLVENT_CONFIGURATION} -p ${TOP}
gmx grompp -f em.mdp -p ${TOP} -c solv.gro -o em
gmx mdrun -v -deffnm em
echo 0 | gmx trjconv -f em.trr -o em.xtc -pbc whole -s em.tpr  # we need to keep all the molecules together without jumping across the box
Last_Frame.py -t em.xtc -g em.gro -o em.gro  # extract the last frame (which is technically the same thing as em.gro output earlier)
gmx grompp -f npt.mdp -p ${TOP} -c em.gro -o ${OUTPUT}
gmx mdrun -v -deffnm ${OUTPUT}
echo "Solvated system set up!"