#!/usr/bin/env bash

# Copy files from a remote resource where there is an equivalent directory structure.

SCP_BASE_PATH="/scratch/summit/beco4952/viscosity_diffusivity"  # past this path, the file structure must be the same as LOCAL_BASE_PATH
LOCAL_BASE_PATH="/home/bcoscia/Documents/Gromacs/viscosity/"  # past this path, the file structure must be the same as SCP_BASE_PATH
SCP_ADDRESS="beco4952@login-new.rc.colorado.edu"
MAX_DEPTH=4

# loop through all subdirectories up to MAX_DEPTH deep
for d in $(find ${PWD}/* -maxdepth ${MAX_DEPTH} -type d); do
  RESOURCE_PATH=${SCP_BASE_PATH}/${d#${LOCAL_BASE_PATH}}  # location of files on HPC resource
  scp ${SCP_ADDRESS}:/${RESOURCE_PATH}/* ${d}  # copy files from HPC resource to local equivalent directory
done