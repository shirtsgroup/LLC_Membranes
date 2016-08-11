#!/usr/bin/env bash

# A script meant for crontab that will backup files on HPC resources by moving them to the correct filesystem

RESOURCE='bridges'

while getopts "r:" opt; do
    case $opt in
    r) RESOURCE=$OPTARG;;
    esac
done

if [ ${RESOURCE} == 'bridges' ]; then
    cp -r /pylon1/ct4s8bp/bjc/* pylon2/ct4s8bp/bjc  # copy files from pylon1 (scratch) to pylon2 (a persistent filesystem)
fi