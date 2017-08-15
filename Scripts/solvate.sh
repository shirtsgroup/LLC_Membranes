#!/usr/bin/env bash

# To get this going, all you need is an initial configuration

WATER_LAYER_THICKNESS=3  # nm of water between membrane layers
GRO="wiggle.gro"  # name of .gro file to solvate
SOLVENT_CONFIGURATION="spc216.gro"
BUILD_MONOMER="NAcarb11V"
TOP="topol.top"
OUTPUT="wiggle"
ref_atoms="C C1 C2 C3 C4 C5"  # atoms used as a reference for calculating thickness based on max/min z positions
monomers_per_layer=5
concentration=0.1 # concentration of salt to add

while getopts "t:g:s:b:p:o:r:m:c:" opt; do
    case $opt in
        t) WATER_LAYER_THICKNESS=$OPTARG;;
        g) GRO=$OPTARG;;
        s) SOLVENT_CONFIGURATION=$OPTARG;;
        b) BUILD_MONOMER=$OPTARG;;
        p) TOP=$OPTARG;;
        o) OUTPUT=$OPTARG;;
        r) ref_atoms=$OPTARG;;
        m) monomer_per_layer=$OPTARG;;
        c) concentration=$OPTARG;;
    esac
done

# put reference atoms in the format of a python list
ref="["
for atom in ${ref_atoms}; do
    ref="${ref}'${atom}',"
done
ref="${ref%?}]"  # %? trims the last character in the string (gets rid of an extra comma)

input.py -b ${BUILD_MONOMER} -c ${GRO} -l 50 --solvate  # create input file

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # the directory where this script is located.
cp ${DIR}/../top/vdwradii.dat .  # copy over the modify vdwradii.dat made for solvation

THICKNESS=$(python -c "from llclib import physical; print(physical.thickness('${GRO}', ${ref})[0])")  # find the membrane thickness

ZBOX=`echo "${THICKNESS} + ${WATER_LAYER_THICKNESS}" | bc -l`  # add the water layer thickness to the membrane thickness

XBOX=$(tail -n 1 ${GRO} | head -c 10)  # get the last line of the gro file to get the box dimensions
# for simplicity, I'm going to assume XBOX only has one component and YBOX=XBOX.

gmx editconf -f ${GRO} -o extended.gro -bt triclinic -c -box ${XBOX} ${XBOX} ${ZBOX} -angles 90 90 120  # create a new box with enough space for water
gmx solvate -cp extended.gro -o solv.gro -cs ${SOLVENT_CONFIGURATION} -p ${TOP}  # solvate box
gmx grompp -f em.mdp -p ${TOP} -c solv.gro -o em

if [ "${concentration}" != 0 ]; then
    # determine the number of ions to add to the system
    conv=0.6022  # avogadros number times 1e-24 (conversion from nm^3 to L)
    V=$(echo "${WATER_LAYER_THICKNESS} * ${XBOX}" | bc -l)  # volume that can be occupied by ions
    n=$(printf %.0f $(echo "${V} * ${conv} * ${concentration}" | bc -l))  # number of ions to add
    conc=$(printf %.3f $(echo "${n} / (${V} * ${conv})" | bc -l))
    echo 'SOL' | gmx genion -s em.tpr -o solv_ions.gro -p ${TOP} -pname NA -nname CL -nn ${n} -np ${n}
    gmx grompp -f em.mdp -p ${TOP} -c solv_ions.gro -o em
fi

gmx mdrun -v -deffnm em  # energy minimize
echo 0 | gmx trjconv -f em.trr -o em.xtc -pbc whole -s em.tpr  # we need to keep all the molecules together without jumping across the box
Last_Frame.py -t em.xtc -g em.gro -o em.gro  # extract the last frame (which is technically the same thing as em.gro output earlier)
gmx grompp -f npt.mdp -p ${TOP} -c em.gro -o ${OUTPUT}
gmx mdrun -v -deffnm ${OUTPUT}

echo "Solvated system set up!"

if [ "${concentration}" != 0 ]; then
    echo "Actual concentration = ${conc} M"
fi