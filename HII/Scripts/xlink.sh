#!/bin/bash
set -e
# A script to iteratively crosslink a LLC system

CUTOFF=2.5  # percent of distribution that will be labeled close enough to bond
TERM_PROB=5  # probability of termination of carbons meeting bonding criteria
GRO='wiggle.gro'  # initial .gro file to be crosslinked
CUTOFF_RAD=20  # percent of distribution that are labeled close enough to bond in the context of reactive radicals
SIM_LENGTH=1  # picoseconds of simulation between iterations
XLINKS=0
FRAMES=50
DEGREE=.9  # Degree of crosslinking
NO_TAILS=3
MONOMER="NAcarb11Vd"
ITERATION=0  # starting iteration
STOP=0
TOP='topol.top'
MDP='npt.mdp'
INIT_CONFIG='initial.gro'
TIMESTEP=0.001
EMSTEPS=5000

while getopts "c:t:g:d:s:x:f:m:C:p:M:i:S:E:T:" opt; do
    case $opt in
    c) CUTOFF=$OPTARG;;
    t) TERM_PROB=$OPTARG;;
    g) GRO=$OPTARG;;
    d) CUTOFF_RAD=$OPTARG;;
    s) SIM_LENGTH=$OPTARG;;
    x) XLINKS=$OPTARG;;
    f) FRAMES=$OPTARG;;
    m) MONOMER=$OPTARG;;
    C) ITERATION=$OPTARG;;  #if the script is being continued
    p) TOP=$OPTARG;;
    M) MDP=$OPTARG;;
    i) INIT_CONFIG=$OPTARG;;
    S) TIMESTEP=$OPTARG;;
    E) EMSTEPS=$OPTARG;;
    T) STOP=$OPTARG;;
    esac
done

export GMX_MAX_BACKUP=-1

if [ ${ITERATION} != 0 ]; then
    ITERATION=$((ITERATION-1))
    XLINKS=$(tail xlink_${ITERATION}.log -n 2 | head -n 1 | cut -c 19-22)
    TERM=$(tail -n 6 xlink_${ITERATION}.log | head -n 1 | cut -c 26-29)
    if [ ${STOP} == 1 ]; then
        :
    else
        STOP=$(tail -n 1 xlink_${ITERATION}.log | cut -c 26)
    fi
    ITERATION=$((ITERATION+1))
fi

input.py -b ${MONOMER} -c ${GRO} -x -l ${SIM_LENGTH} -f ${FRAMES} --genvel yes -d ${TIMESTEP} -s ${EMSTEPS}

while [ ${STOP} -eq 0 ]; do
    if [ ${ITERATION} == 0 ]; then
        xlink.py -i ${GRO} -c ${CUTOFF}  -e ${TERM_PROB} -r ${ITERATION} -d ${CUTOFF_RAD} -x ${XLINKS} -m ${MONOMER} --nogap
    else
        xlink.py -i ${GRO} -c ${CUTOFF}  -e ${TERM_PROB} -r ${ITERATION} -d ${CUTOFF_RAD} -y crosslinked_new.itp -x ${XLINKS} --nogap
    fi

    cp crosslinked_new.itp crosslinked_${ITERATION}.itp
    mv xlink.log xlink_${ITERATION}.log

    if [ ${ITERATION} == 0 ]; then
        gmx grompp -f em.mdp -p topol.top -c ${GRO} -o em
    else
        gmx grompp -f em.mdp -p topol.top -c wiggle.gro -o em
    fi
    gmx mdrun -v -deffnm em
    gmx grompp -f ${MDP} -p topol.top -c em.gro -o wiggle
    gmx mdrun -v -deffnm wiggle
#    echo 0 | gmx trjconv -f wiggle.gro -s wiggle.gro -pbc atom -o wiggle.gro -ur tric
    XLINKS=$(tail xlink_${ITERATION}.log -n 2 | head -n 1 | cut -c 19-22)
    TERM=$(tail -n 6 xlink_${ITERATION}.log | head -n 1 | cut -c 26-29)
    STOP=$(tail -n 1 xlink_${ITERATION}.log | cut -c 26)
    ITERATION=$((ITERATION+1))
    echo ${TERM}
#    files=$(ls ./\#*.cache 2> /dev/null | wc -l)
#    if [ **"$files" != "0"** ]; then
#        rm \#*
#    fi
done

xlink.py -i ${GRO} -c ${CUTOFF}  -e ${TERM_PROB} -r ${ITERATION} -d ${CUTOFF_RAD} -y crosslinked_new.itp -x ${XLINKS} -S 'yes' --nogap
# Put all the atoms back in the unit cell
#gmx trjconv -f wiggle.gro -o unitcell.gro -ur tric -pbc atom -s ${INIT_CONFIG}  # good for visualization, bad for simulation (inconsistent shifts)

cp wiggle.gro wiggle_dummies.gro

Cleanup_Top.py  # get rid of dummy atoms in .itp and .gro

find ${TOP} -type f -exec sed -i 's/crosslinked_new.itp/crosslinked.itp/g' {} \;

gmx grompp -f em.mdp -p ${TOP} -c wiggle_no_dummies.gro -o em
gmx mdrun -v -deffnm em
gmx grompp -f npt.mdp -p ${TOP} -c em.gro -o wiggle
gmx mdrun -v -deffnm wiggle

mail -s "Crosslinking Done" -a xlink.log benjamin.coscia@colorado.edu <<< "The crosslinking has terminated, woo!"

rm \#*  # get rid of backup files
