#!/bin/bash

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
B_ALLOC='ct4s8bp'  # bridges allocation
GMX_LOC="/home/bjc/Programs/Gromacs/bin/GMXRC"  # location of gmx executable if you have compiled your own

while getopts "e:n:t:c:l:f:j:x:r:h:m:s:a:g" opt; do
    case ${opt} in
    e) EXTENSION=$OPTARG;;
    n) NODES=$OPTARG;;
    t) TPR=$OPTARG;;
    c) CPT=$OPTARG;;
    l) LOG=$OPTARG;;
    j) JOB_ARRAY+=($OPTARG);;  # optional since it can be made automatically
    x) EXCLUSIONS+=($OPTARG);;  # optional. Use to exclude jobs from JOB_ARRAY the program automatically produces JOB_ARRAY
    # Note: to have multiple exlusions, use the -x flags multiple times -x job1 -x job2 -x job3 etc.
    r) RESOURCE=$OPTARG;;
    h) HOURS=$OPTARG;;
    m) MIN=$OPTARG;;
    s) SEC=$OPTARG;;
    a) B_ALLOC=$OPTARG;;
    g) GMX_LOC=$OPTARG;;
    esac
done

# check what machine this script is being run on so we get the username correct, request the correct partition and make
# sure that we are allowed to run for the amount of hours requested

if [ ${RESOURCE} == 'janus' ]; then
    USER='beco4952'  # Can't use environment variable because of cron
    NTASKSPERNODE=1
    NP=$((NODES*2))
    OUTPUT="/lustre/janus_scratch/${USER}/Cron_Output"  # place where files like job_array, queue and dir will be stored
    GMX_LOC="/projects/${USER}/bin/GMXRC"
    if (( ${HOURS}>24 )); then
        QOS='janus-long'
    else
        QOS='janus'
    fi
    if (( ${HOURS}>168 )); then
        echo 'That is too long'
        exit  # Stop the script
    fi
fi

if [ ${RESOURCE} == 'bridges' ]; then
    USER='bjc'
    QOS='RM'
    NTASKSPERNODE=4
    NP=$((NTASKSPERNODE*NODES))
    OUTPUT="/pylon1/${B_ALLOC}/${USER}/Cron_Output"
    GMX_LOC="/home/${USER}/Programs/Gromacs/bin/GMXRC"
    if (( ${HOURS}>48 )); then
        echo 'That is too long'
        exit
    fi
fi

# First check that the job is running. It may be stuck in the queue
# Create an array of jobs:

squeue --user=$USER > ${OUTPUT}/queue  # save the queue output to a file
NO_LINES=$(wc -l < ${OUTPUT}/queue)  # find the number of lines in queue
for i in $(seq 2 $NO_LINES); do  # add the job id for all jobs in the user's queue to JOB_ARRAY. Exclude the first line
    echo $(sed "${i}q;d" ${OUTPUT}/queue | cut -c 12-18) >> ${OUTPUT}/job_array
done

sort ${OUTPUT}/job_array | uniq > ${OUTPUT}/job_array_uniq # get rid of duplicates

JOB_ARRAY=()  # intialize job array
ENTRIES=$(wc -l < ${OUTPUT}/job_array_uniq)  # find the number of lines in queue
for i in $(seq 1 ${ENTRIES}); do  # add the job id for all jobs in the user's queue to JOB_ARRAY. Exclude the first line
    JOB_ARRAY+=($(sed "${i}q;d" ${OUTPUT}/job_array_uniq | cut -c 1-7))
done

# Remove any exclusions
for i in ${EXCLUSIONS[@]}; do
   JOB_ARRAY=( "${JOB_ARRAY[@]/$i}" )
done

for JOB_ID in ${JOB_ARRAY[@]}; do  # look at all running jobs
    squeue --user=$USER > ${OUTPUT}/queue  # save the queue information to a temporary file called queue
    LINE=$(sed -n "/${JOB_ID}/=" ${OUTPUT}/queue) # find the line containing the job_id of interest

    if [[ -z ${LINE} ]]; then  # if $LINE is an empty variable, meaning it doesn't exist in the queue, then the job is probably done
        DONE=1
    else
        STATUS=$(sed "${LINE}q;d" ${OUTPUT}/queue | cut -b 48,49)  # search that line to find the job's status

        if [ ${STATUS} != ' R' ]; then  # if the job is not running but in the queue or completing, leave it. Exit the whole script
            exit
        fi
    fi

    # This point in the loop is reached only if the job is running or finished

    # figure out the working directory for the job ID
    scontrol show job ${JOB_ID} > ${OUTPUT}/dir  # save all the information pertaining to the job id to a file
    WORKDIRLINE=$(sed -n '/WorkDir/=' ${OUTPUT}/dir)  # find the line containing the working directory for the job
    WORKDIRFULL=$(sed "${WORKDIRLINE}q;d" ${OUTPUT}/dir)  # save the entire line
    WORKDIR=$(sed "${WORKDIRLINE}q;d" ${OUTPUT}/dir | cut -c 12-${#WORKDIRFULL})  # cut from the 12th character to the length of the line

    # We need to monitor the currently running job to see if it is still going, or if it is stalled/finished
    tail ${WORKDIR}/${LOG} > ${WORKDIR}/log_${JOB_ID}  # Save the last 10 lines of the .log file to log_prev. This will be used as a 1st comparison

    CHECK=0

    # We need a log file from a previous iteration for comparison. If the job just started, then one does not
    # exist. I make a dummy file below so that log_${JOB_ID} and log_${JOB_ID}_prev will be different which will prevent
    # the script from extending and restarting the job every time
    if [ -e "${WORKDIR}/log_${JOB_ID}_prev" ]; then  # check for the existence of the file log_${JOB_ID}_prev
        :  # if it exists, do nothing
    else
        echo 'Placeholder' > ${WORKDIR}/log_${JOB_ID}_prev  # if it doesn't exist, make a place holding log_${JOB_ID}_prev file
    fi

    cmp ${WORKDIR}/log_${JOB_ID} ${WORKDIR}/log_${JOB_ID}_prev || export CHECK=1  # if the two files are different then CHECK=1 ... this is what we want unless the job
    # is finished. Otherwise (CHECK=0), it is an indication of a job that needs to be restarted (a stalled or finished job)
    # now change log to the new log_prev for the next iteration
    mv ${WORKDIR}/log_${JOB_ID} ${WORKDIR}/log_${JOB_ID}_prev

    echo ${CHECK}

    if [ $CHECK == 0 ]; then  # if the files are the same (meaning output is no longer being written)

        if [ $DONE == 1 ]; then
            :  # indicates a job that has completed (i.e. no longer in the queue)
        else
            scancel ${JOB_ID}  # run this command to cancel the job only if the job is still running
        fi

        if [ ${RESOURCE} == 'janus' ]; then  # write a janus batch submission script
            echo "#! /bin/bash" > ${WORKDIR}/Extend_Sim_new.sh
            echo >> ${WORKDIR}/Extend_Sim_new.sh
            echo "#SBATCH --qos ${QOS}" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "#SBATCH --nodes ${NODES}" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "#SBATCH --ntasks-per-node ${NTASKSPERNODE}" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "#SBATCH --time ${HOURS}:${MIN}:${SEC}" >> ${WORKDIR}/Extend_Sim_new.sh
            echo >> ${WORKDIR}/Extend_Sim_new.sh
            echo "ml slurm" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "ml gromacs" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "ml python/2.7.10" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "ml numpy" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "source ${GMX_LOC}" >> ${WORKDIR}/Extend_sim_new.sh
            echo >> ${WORKDIR}/Extend_Sim_new.sh
            echo "gmx convert-tpr -s ${TPR} -extend ${EXTENSION} -o ${TPR}" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "mpirun -np ${NP} gmx_mpi mdrun -s ${TPR} -cpi ${CPT} -append -v -deffnm "${TPR//.tpr}"" >> ${WORKDIR}/Extend_Sim_new.sh
            sed -i "/${JOB_ID}/d" -i ${OUTPUT}/job_array # remove all instances of this job_id from job_array
        fi

        if [ ${RESOURCE} == 'bridges' ]; then  # write a bridges batch submission script
            echo "#!/bin/bash" > ${WORKDIR}/Extend_Sim_new.sh
            echo >> ${WORKDIR}/Extend_Sim_new.sh
            echo "#SBATCH --partition=${QOS}" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "#SBATCH -N ${NODES} --tasks-per-node=${NTASKSPERNODE}" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "#SBATCH -t ${HOURS}:${MIN}:${SEC}" >> ${WORKDIR}/Extend_Sim_new.sh
            echo >> ${WORKDIR}/Extend_Sim_new.sh
            # This block is optional
            echo "set echo" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "set -x" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "source /usr/share/Modules/init/bash" >> ${WORKDIR}/Extend_Sim_new.sh
            echo >> ${WORKDIR}/Extend_Sim_new.sh
            # End optional part
            # Add other modules if you need them
            echo "module load gromacs" >> ${WORKDIR}/Extend_Sim_new.sh
            echo source ${GMX_LOC} >> ${WORKDIR}/Extend_Sim_new.sh
            echo >> ${WORKDIR}/Extend_Sim_new.sh
            echo "cd $SLURM_SUBMIT_DIR" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "echo "$SLURM_NPROCS=" $SLURM_NPROCS" >> ${WORKDIR}/Extend_Sim_new.sh
            echo >> ${WORKDIR}/Extend_Sim_new.sh
            echo "gmx convert-tpr -s ${TPR} -extend ${EXTENSION} -o ${TPR}" >> ${WORKDIR}/Extend_Sim_new.sh
            echo "mpirun -np ${NP} gmx_mpi mdrun -s ${TPR} -cpi ${CPT} -append -v -deffnm "${TPR//.tpr}"" >> ${WORKDIR}/Extend_Sim_new.sh
            sed -i "/${JOB_ID}/d" -i ${OUTPUT}/job_array # remove all instances of this job_id from job_array
        fi
    fi
done
