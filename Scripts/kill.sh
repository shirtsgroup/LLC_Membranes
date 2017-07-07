#!/bin/bash

squeue --user=$USER > jobs # find out what jobs are running

jobs=$(get_job_ids.py -i jobs)  # get the job ids

for i in $jobs; do  # check all jobs
	scontrol show job $i > job_info; # get info about running jobs
	path=$(get_job_path.py -i job_info);  # get paths to the StdOut
	if [ -f ${path} ]; then  # check if the file at the given path exists. If it doesn't exist then job hasn't started running yet
		if grep -q 'LINCS' ${path}; then # Check if there are LINCS warnings
			scancel $i  # cancel the job because it is probably stopped
		fi
	fi
done
