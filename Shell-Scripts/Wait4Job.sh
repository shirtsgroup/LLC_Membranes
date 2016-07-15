#!/bin/bash

WORKINGDIR=$PWD  # this will probably always be the case
FILE='Extend_Sim.sh'  # this will also probably always be the case
SLEEP=180  # default sleep for 3 minutes

while getopts "p:f:s:" opt; do
    case ${opt} in
    p) PATH=$OPTARG;;
    f) FILE=$OPTARG;;
    s) SLEEP=$OPTARG;;
    esac
done

while :  # starts an infinite loop
do
	until [ -f ${WORKINGDIR}/${FILE} ]  # runs this loop until ${FILE} exists in ${WORKINGDIR}
	do
		sleep ${SLEEP}  # sleep for a specified amount of time between checks
	done
	sbatch ${WORKINGDIR}/Extend_Sim.sh  # submit the job
done  # loop never ends so that is continuously checks for new Extend_Sim_new.sh files