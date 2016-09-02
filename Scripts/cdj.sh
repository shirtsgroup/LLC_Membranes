#!/bin/bash

# Changes directories to where a job associated with the input job ID is running
# Run this script prefixed by '.' so that it actually changes directories, otherwise directories are changed in a bash
# subshell and the user will see no change

while getopts "j:" opt; do
    case $opt in
    j) JOBID=$OPTARG;;
    esac
done

scontrol show job ${JOBID} > job

workdir=$(sed -n 's/^   WorkDir=//p' job)

cd $workdir
