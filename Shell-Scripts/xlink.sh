#!/bin/bash
set -e
# A script to iteratively crosslink a LLC system

N_ITER=10  # number of crosslinking iterations
CUTOFF=5  # percent of distribution that will be labeled close enough to bond
TERM_PROB=5  # probability of termination of carbons meeting bonding criteria
GRO='wiggle.gro'  # initial .gro file to be crosslinked
CUTOFF_RAD=10  # percent of distribution that are labeled close enough to bond in the context of reactive radicals
SIM_LENGTH=.01  # picoseconds of simulation between iterations
XLINKS=0

while getopts "n:c:t:g:d:s:x:" opt; do
    case $opt in
    n) N_ITER=$OPTARG;;
    c) CUTOFF=$OPTARG;;
    t) TERM_PROB=$OPTARG;;
    g) GRO=$OPTARG;;
    d) CUTOFF_RAD=$OPTARG;;
    s) SIM_LENGTH=$OPTARG;;
    x) XLINKS=$OPTARG;;
    esac
done

ITERATION=0  # starting iteration
TERM=0  # starting number of vinyl groups that have been terminated (both c1 and c2)

Write_Input.py -x on -L ${SIM_LENGTH} -D 0.001 # -I cg

while [  ${TERM} -lt 1440 ]; do
    if [ ${ITERATION} == 0 ]; then
        xlink.py -i ${GRO} -c ${CUTOFF}  -e ${TERM_PROB} -r ${ITERATION} -d ${CUTOFF_RAD} -x ${XLINKS}
    else
        xlink.py -i ${GRO} -c ${CUTOFF}  -e ${TERM_PROB} -r ${ITERATION} -d ${CUTOFF_RAD} -y crosslinked_new.itp -x ${XLINKS}
    fi
    cp crosslinked_new.itp crosslinked_${ITERATION}.itp
    mv xlink.log xlink_${ITERATION}.log
    gmx grompp -f em.mdp -p NaPore.top -c wiggle.gro -o em
    gmx mdrun -v -deffnm em
    gmx grompp -f wiggle.mdp -p NaPore.top -c em.gro -o wiggle
    gmx mdrun -v -deffnm wiggle
    XLINKS=$(tail xlink_${ITERATION}.log -n 1 | cut -c 19-22)
    TERM=$(tail -n 5 xlink_9.log | head -n 1 | cut -c 26-29)
    ITERATION=$((ITERATION+1))
    echo ${TERM}
done

rm \#*