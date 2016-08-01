#!/bin/bash
set -e
# A script to iteratively crosslink a LLC system

N_ITER=10  # number of crosslinking iterations
CUTOFF=5  # percent of distribution that will be labeled close enough to bond
TERM_PROB=5  # probability of termination of carbons meeting bonding criteria
GRO='wiggle.gro'  # initial .gro file to be crosslinked
CUTOFF_RAD=10  # percent of distribution that are labeled close enough to bond in the context of reactive radicals
SIM_LENGTH=.01  # picoseconds of simulation between iterations

while getopts "n:c:t:g:d:s:" opt; do
    case $opt in
    n) N_ITER=$OPTARG;;
    c) CUTOFF=$OPTARG;;
    t) TERM_PROB=$OPTARG;;
    g) GRO=$OPTARG;;
    d) CUTOFF_RAD=$OPTARG;;
    s) SIM_LENGTH=$OPTARG;;
    esac
done

ITERATION=0

Write_Input.py -x on -L ${SIM_LENGTH} -D 0.001 # -I cg

for i in $(seq 0 $((N_ITER-1))); do
    if [ ${ITERATION} == 0 ]; then
        xlink.py -i ${GRO} -c ${CUTOFF}  -e ${TERM_PROB} -r ${ITERATION} -d ${CUTOFF_RAD}
    else
        xlink.py -i ${GRO} -c ${CUTOFF}  -e ${TERM_PROB} -r ${ITERATION} -d ${CUTOFF_RAD} -y crosslinked_new.itp
    fi
    cp crosslinked_new.itp crosslinked_${ITERATION}.itp
    mv xlink.log xlink_${ITERATION}.log
    gmx grompp -f em.mdp -p NaPore.top -c wiggle.gro -o em
    gmx mdrun -v -deffnm em
    gmx grompp -f wiggle.mdp -p NaPore.top -c em.gro -o wiggle
    gmx mdrun -v -deffnm wiggle
    ITERATION=$((ITERATION+1))
done

rm \#*