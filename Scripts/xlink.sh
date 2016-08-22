#!/bin/bash
set -e
# A script to iteratively crosslink a LLC system

CUTOFF=5  # percent of distribution that will be labeled close enough to bond
TERM_PROB=5  # probability of termination of carbons meeting bonding criteria
GRO='wiggle.gro'  # initial .gro file to be crosslinked
CUTOFF_RAD=10  # percent of distribution that are labeled close enough to bond in the context of reactive radicals
SIM_LENGTH=.01  # picoseconds of simulation between iterations
XLINKS=0
FRAMES=50
DEGREE=.9  # Degree of crosslinking
NO_MONOMERS=480
NO_TAILS=3
MONOMER="monomer2"
ITERATION=0  # starting iteration
STOP=0

while getopts "c:t:g:d:s:x:f:m:C:" opt; do
    case $opt in
    c) CUTOFF=$OPTARG;;
    t) TERM_PROB=$OPTARG;;
    g) GRO=$OPTARG;;
    d) CUTOFF_RAD=$OPTARG;;
    s) SIM_LENGTH=$OPTARG;;
    x) XLINKS=$OPTARG;;
    f) FRAMES=$OPTARG;;
    m) MONOMER=$OPTARG;;
    C) ITERATION=$OPTARG;;  #if the sript is being continued
    esac
done

if [ ${ITERATION} != 0 ]; then
    ITERATION=$((ITERATION-1))
    XLINKS=$(tail xlink_${ITERATION}.log -n 2 | head -n 1 | cut -c 19-22)
    TERM=$(tail -n 6 xlink_${ITERATION}.log | head -n 1 | cut -c 26-29)
    STOP=$(tail -n 1 xlink_${ITERATION}.log | cut -c 26)
    ITERATION=$((ITERATION+1))
fi

Write_Input.py -x on -L ${SIM_LENGTH} -D 0.001 -f ${FRAMES} # -I cg

while [ ${STOP} -eq 0 ]; do
    if [ ${ITERATION} == 0 ]; then
        xlink.py -i ${GRO} -c ${CUTOFF}  -e ${TERM_PROB} -r ${ITERATION} -d ${CUTOFF_RAD} -x ${XLINKS} -m ${MONOMER}
    else
        xlink.py -i ${GRO} -c ${CUTOFF}  -e ${TERM_PROB} -r ${ITERATION} -d ${CUTOFF_RAD} -y crosslinked_new.itp -x ${XLINKS}
    fi
    cp crosslinked_new.itp crosslinked_${ITERATION}.itp
    mv xlink.log xlink_${ITERATION}.log
    gmx grompp -f em.mdp -p NaPore.top -c wiggle.gro -o em
    gmx mdrun -v -deffnm em
    gmx grompp -f wiggle.mdp -p NaPore.top -c em.gro -o wiggle
    gmx mdrun -v -deffnm wiggle
    XLINKS=$(tail xlink_${ITERATION}.log -n 2 | head -n 1 | cut -c 19-22)
    TERM=$(tail -n 6 xlink_${ITERATION}.log | head -n 1 | cut -c 26-29)
    STOP=$(tail -n 1 xlink_${ITERATION}.log | cut -c 26)
    ITERATION=$((ITERATION+1))
    echo ${TERM}
    if test -f "\#*"; then rm \#*; fi  # remove backups if they exist
done

xlink.py -i ${GRO} -c ${CUTOFF}  -e ${TERM_PROB} -r ${ITERATION} -d ${CUTOFF_RAD} -y crosslinked_new.itp -x ${XLINKS} -S 'yes'

Cleanup_Top.py  # get rid of dummy atoms in .itp and .gro

find NaPore.top -type f -exec sed -i 's/crosslinked_new.itp/crosslinked.itp/g' {} \;

gmx grompp -f em.mdp -p NaPore.top -c wiggle_no_dummies.gro -o em
gmx mdrun -v -deffnm em
gmx grompp -f wiggle.mdp -p NaPore.top -c em.gro -o wiggle
gmx mdrun -v -deffnm wiggle

mail -s "Crosslinking Done" -a xlink.log benjamin.coscia@colorado.edu <<< "The crosslinking has terminated, woo!"

rm \#*  # get rid of backup files
