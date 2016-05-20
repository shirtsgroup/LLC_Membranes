#!/bin/bash

# Script to copy trajectory file from Janus, convert it's trajectory and do analysis

FOLDER='Output'
PATH='/projects/beco4952/Gromacs/Pores/
# Flags
# -f  : FOLDER ... folder to be copied over (should include at least .trr, .tpr files)
# -p  : PATH ... path to file/folder which will be copied over

while getopts "f:p:" opt; do
    case $opt in
    f) FOLDER=$OPTARG;;
    p) PATH=$OPTARG;;
    esac
done

scp -r $FOLDER beco4952@login.rc.colorado.edu:$PATH $PWD  # copy folder to directory this script is run in

cd $FOLDER

gmx trjconv -f wiggle.trr -s wiggle.tpr -o wiggle_traj.gro

Cylindricity_Traj.py -i wiggle_traj.gro

