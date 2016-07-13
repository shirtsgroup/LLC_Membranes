#!/bin/bash

WORKINGDIR=$PWD
FILE='Extend_Sim_new.sh'
while getopts "p:f:" opt; do
    case ${opt} in
    p) PATH=$OPTARG;;
    f) FILE=$OPTARG;;
    esac
done

while :
do
	until [ -f ${WORKINGDIR}/${FILE} ]
	do
		sleep 5
	done
	mv ${WORKINGDIR}/${FILE} ${WORKINGDIR}/Extend_Sim.sh
	sbatch Extend_Sim.sh
done
