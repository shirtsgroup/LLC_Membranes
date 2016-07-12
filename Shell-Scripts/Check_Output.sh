#!/usr/bin/env bash

EXTENSION=50  # number of picoseconds that you want to extend the simulation by
TPR="wiggle.tpr"  # .tpr file to be edited for simulation extension
CPT="wiggle.cpt"  # .cpt file to be referenced for simulation extension
LOG="wiggle.log"  # .log file to be monitored for stalling/ending of simulations
RESOURCE="janus"  # name of machine -- all lowercase
EXCLUSIONS=()  # job id's to be excluded (i.e. other jobs that are running). Make this an environment variable
HOURS=48  # hours requested for sim extension. Default 48 since that's the min of the max queues on janus & bridges
MIN=00  # minutes requested for sim extension
SEC=00  # seconds requested for sim extension
NODES=8  # number of nodes requested

while getopts "e:n:t:c:l:f:j:x:r:h:m:s" opt; do
    case $opt in
    e) EXTENSION=$OPTARG;;
    t) TPR=$OPTARG;;
    c) CPT=$OPTARG;;
    l) LOG=$OPTARG;;
    j) JOB_ARRAY+=($OPTARG);;  # optional since it can be made automatically
    x) EXCLUSIONS+=($OPTARG);;  # optional. Use to exclude jobs from JOB_ARRAY the program automatically produces JOB_ARRAY
    r) RESOURCE=$OPTARG;;
    h) HOURS=$OPTARG;;
    m) MIN=$OPTARG;;
    s) SEC=$OPTARG;;
    esac
done

# check what machine this script is being run on so we get the username correct, request the correct partition and make
# sure that we are allowed to run for the amount of hours requested

if [ ${RESOURCE}=='janus' ]; then
    USER='beco4952'
    NTASKSPERNODE=1
    NP=$((NODES*2))
    if [ ${HOURS}>24 ]; then
        QOS='janus-long'
    else
        QOS='janus'
    fi
    if [ ${HOURS}>168 ]; then
        echo 'That is too long'
        exit
    fi
elif [ ${RESOURCE}=='bridges' ]; then
    USER='bjc'
    QOS='RM'
    NTASKSPERNODE=28
    NP=$((NTASKSPERNODE*NODES))
    if [ ${HOURS}>48 ]; then
        echo 'That is too long'
        exit
    fi
fi

# First check that the job is running. It may be stuck in the queue
# Submit jobs and simultaneously add them to the JOB_ARRAY:
# JOB_ARRAY+=($(sbatch Run_Janus.sh | cut -b 21,22,23,24,25,26,27))
# or create an array automatically:

JOB_ARRAY=()  # initialize JOB_ARRAY
squeue --user=$USER > queue  # save the queue output to a file
NO_LINES=$(wc -l < queue)  # find the number of lines in queue
for i in ${seq 2 $NO_LINES}; do  # add the job id for all jobs in the user's queue to JOB_ARRAY. Exclude the first line
    JOB_ARRAY+=($(sed "${i}q;d" queue | cut -c 12-18))

# Remove any exclusions
for i in ${EXCLUSIONS[@]}; do
   JOB_ARRAY=( "${JOB_ARRAY[@]/$i}" )
done

ITERATION=1
for JOB_ID in ${JOB_ARRAY[@]}; do  # look at all running jobs
    squeue --user=$USER > queue  # save the queue information to a temporary file called queue
    LINE=$(sed -n "/${JOB_ID}/=" queue)  # find the line containing the job_id of interest
    if [[ -z ${LINE} ]]; then  # if $LINE is empty, meaning it doesn't exist in the queue, then the job is probably done
        : # do nothing and move on to extending the job
    else
        STATUS=$(sed "${LINE}q;d" queue | cut -b 48,49)  # search that line to find the job's status

        if [ ${STATUS} != ' R' ]; then  # if the job is not running but in the queue or completing, leave it. Exit the whole script
            exit
        fi
    fi

    # This point in the loop is reached only if the job is running or finished

    # figure out the working directory for the job ID
    scontrol show job ${JOB_ID} > dir  # save all the information pertaining to the job id to a file
    WORKDIRLINE=$(sed -n '/WorkDir/=' dir)  # find the line containing the working directory for the job
    WORKDIRFULL=$(sed "${WORKDIRLINE}q;d" dir)  # save the entire line
    WORKDIR=$(sed "${WORKDIRLINE}q;d" dir | cut -c 12-${#WORKDIRFULL})  # cut from the 12th character to the length of the line

    # We need to monitor the currently running job to see if it is still going, or if it is stalled/finished
    tail ${WORKDIR}/${LOG} > log_${JOB_ID}  # Save the last 10 lines of the .log file to log_prev. This will be used as a 1st comparison

    CHECK=0

    # We need a log file from a previous iteration for comparison. If this is the first iteration, then one does not
    # exist. I make a dummy file below so that log_${JOB_ID} and log_${JOB_ID}_prev will be different which will prevent
    # the script from extending and restarting the job every time
    if [ ${ITERATION}==1 ]; then  # if it's the first iteration
        echo 'Placeholder' > log_${JOB_ID}_prev  # make a placeholding log_${JOB_ID}_prev
    fi

    cmp log_${JOB_ID} log_${JOB_ID}_prev || export CHECK=1  # if the two files are different, CHECK=1 ... this is what we want unless the job
    # is finished. Otherwise (CHECK=0), it is an indication of a job that needs to be restarted (a stalled or finished job)
    # now change log to the new log_prev for the next iteration
    mv log_${JOB_ID} log_${JOB_ID}_prev

    if $CHECK == 0; then  # if the files are the same (meaning output is no longer being written)

        scancel ${JOB_ID}  # cancel the job which will be extended
        gmx_mpi convert-tpr -s ${TPR} -extend ${EXTENSION} -o ${TPR}  # make a new, extended .tpr file

        if [ ${RESOURCE}=='janus' ]; then  # write a janus batch submission script
            echo "#! /bin/bash" > ${WORKDIR}/Extend_Sim.sh
            echo >> ${WORKDIR}/Extend_Sim.sh
            echo "#SBATCH --qos ${QOS}" >> ${WORKDIR}/Extend_Sim.sh
            echo "#SBATCH --nodes ${NODES}" >> ${WORKDIR}/Extend_Sim.sh
            echo "#SBATCH --ntasks-per-node ${NTASKSPERNODE}" >> ${WORKDIR}/Extend_Sim.sh
            echo "#SBATCH --time ${HOURS}:${MIN}:${SEC}" >> ${WORKDIR}/Extend_Sim.sh
            echo >> ${WORKDIR}/Extend_Sim.sh
            echo "ml slurm" >> ${WORKDIR}/Extend_Sim.sh
            echo "ml gromacs" >> ${WORKDIR}/Extend_Sim.sh
            echo "ml python/2.7.10" >> ${WORKDIR}/Extend_Sim.sh
            echo "ml numpy" >> ${WORKDIR}/Extend_Sim.sh
            echo >> ${WORKDIR}/Extend_Sim.sh
            echo "mpirun -np ${NP} gmx_mpi mdrun -s ${TPR} -cpi ${CPT} -v -deffnm "${TPR//.tpr}"" >> ${WORKDIR}/Extend_Sim.sh
            sbatch Extend_Sim.sh  # submit job
        elif [ ${RESOURCE}=='bridges' ]; then  # write a bridges batch submission script
            echo "#!/bin/bash" > ${WORKDIR}/Extend_Sim.sh
            echo >> ${WORKDIR}/Extend_Sim.sh
            echo "#SBATCH --partition=${QOS}" >> ${WORKDIR}/Extend_Sim.sh
            echo "#SBATCH -N ${NODES} --tasks-per-node=${NTASKSPERNODE}" >> ${WORKDIR}/Extend_Sim.sh
            echo "#SBATCH -t ${HOURS}:${MIN}:${SEC}" >> ${WORKDIR}/Extend_Sim.sh
            echo >> ${WORKDIR}/Extend_Sim.sh
            echo "set echo" >> ${WORKDIR}/Extend_Sim.sh
            echo "set -x" >> ${WORKDIR}/Extend_Sim.sh
            echo "source /usr/share/Modules/init/bash" >> ${WORKDIR}/Extend_Sim.sh
            echo >> ${WORKDIR}/Extend_Sim.sh
            echo "module load gromacs" >> ${WORKDIR}/Extend_Sim.sh
            echo >> ${WORKDIR}/Extend_Sim.sh
            echo "cd $SLURM_SUBMIT_DIR" >> ${WORKDIR}/Extend_Sim.sh
            echo "echo "$SLURM_NPROCS=" $SLURM_NPROCS" >> ${WORKDIR}/Extend_Sim.sh
            echo >> ${WORKDIR}/Extend_Sim.sh
            echo "mpirun -np ${NP} gmx_mpi mdrun -s ${TPR} -cpi ${CPT} -v -deffnm "${TPR//.tpr}"" >> ${WORKDIR}/Extend_Sim.sh
            sbatch Extend_Sim.sh # submit job
        fi
    fi
    ITERATION=$((ITERATION+1))  # increment iteration variable
done