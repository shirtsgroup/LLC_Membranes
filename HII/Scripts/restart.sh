#!/usr/bin/env bash

squeue --user=$USER > jobs # find out what jobs are running

jobs=$(get_job_ids.py -i jobs)  # get the job ids

for i in $jobs; do
    echo $i >> jobs_running  # list of current jobs running when this script was run
done

# make sure there is an old list of jobs to compare to. This is necessary so it works on the first run
if [ -f jobs_running_old ]; then
    restarts=$(restarts.py)
else
    exit
fi

for i in $restarts; do  # check all jobs that stopped running since last checked
	scontrol show job $i > job_info; # get info about running jobs
	path=$(get_job_path.py -i job_info);  # get paths to the StdOut
	if [ -f ${path} ]; then  # check if the file at the given path exists. If it doesn't exist then job hasn't started running yet
		if grep -q 'LINCS' ${path}; then # Check if there are LINCS warnings
			DIR=$(dirname "${path}") # directory where job is located
			KEEP=$(grep -n 'module' ${DIR}/Run.sh | tail -n1 | cut -d: -f1)  # keep up to this line from the original Run.sh file
			head -n ${KEEP} ${DIR}/Run.sh > ${DIR}/restart.sh  # copy lines up to KEEP into a new batch script
			processes=$(get_node_config.py -i ${DIR}/Run.sh)  # number of processes
			echo '' >> ${DIR}/restart.sh
			echo "mpirun -np ${processes} gmx_mpi mdrun -s wiggle.tpr -cpi wiggle.cpt -append -v -deffnm wiggle" >> ${DIR}/continue.sh # command to continue where the simulation left off
			sbatch --workdir ${DIR} ${DIR}/restart.sh  # submit job for continuation
			time=$(date '+%H:%M')
			mail -s "Job ${i} restarted at $time" benjamin.coscia@colorado.edu <<< "Path to job failure: ${path}"  # email myself the location of failed job
		fi
	fi
done

mv jobs_running jobs_running_old