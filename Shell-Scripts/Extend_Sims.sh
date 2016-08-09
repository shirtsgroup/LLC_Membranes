#!/usr/bin/env bash

# String together jobs on janus and bridges

EXTENSION=50  # number of picoseconds that you want to extend the simulation by
N_EXT=1  # number of extensions
TPR="wiggle.tpr"  # .tpr file to be edited for simulation extension
CPT="wiggle.cpt"  # .cpt file to be referenced for simulation extension
LOG="wiggle.log"  # .log file to be monitored for stalling/ending of simulations
FREQ_M=5  # frequency for checking log file (minutes)

while getopts "e:n:t:c:l:f" opt; do
    case $opt in
    e) EXTENSION=$OPTARG;;
    n) N_EXT=$OPTARG;;
    t) TPR=$OPTARG;;
    c) CPT=$OPTARG;;
    l) LOG=$OPTARG;;
    f) FREQ_M=$OPTARG;;
    esac
done

for ((i=0;i<${N_EXT};i++)); do  # run the simulation and extend it a specified number of times

    # We need to monitor the currently running job to see if it is still going, or if it is stalled/finished
    tail ${LOG} > log_prev  # Save the last 10 lines of the .log file to log_prev. This will be used as a 1st comparison

    crontab -l | { cat; echo "*/${FREQ_M} * * * * tail ${PWD}/${LOG} > ${PWD}/log"  # set up a cron job that will save
    # the last ten lines of the log file as it exists at that time every FREQ_M number of minutes

    export CHECK=0  # if CHECK=0 after the following loop, that means the job is done or stalled


    gmx convert-tpr -s ${TPR} -extend ${EXTENSION} -o ${TPR}

    gmx mdrun -s ${TPR} -cpi ${CPT} -v -deffnm "${TPR//.tpr}"

done