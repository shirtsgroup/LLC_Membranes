#!/usr/bin/env bash

# Script to copy trajectory file from Janus, convert it's trajectory and do analysis

FILE1='wiggle.tpr'
FILE2='wiggle.trr'
PATH2FILE='/projects/beco4952/Gromacs/Pores'

# Flags
# -f  : FILE1 ... .tpr file
# -F  : FILE2 ... .trr file
# -p  : PATH2FILE ... path to file/folder which will be copied over

while getopts "f:F:p:" opt; do
    case $opt in
    f) FILE1=$OPTARG;;
    F) FILE2=$OPTARG;;
    p) PATH2FILE=$OPTARG;;
    esac
done

scp beco4952@login.rc.colorado.edu:$PATH2FILE/\{$FILE1,$FILE2\} ./
#scp beco4952@login.rc.colorado.edu:$PATH/\{$FILE1,$FILE2\} ./ # copy files to directory this script is run in

echo "0" | gmx trjconv -f $FILE2 -s $FILE1 -o wiggle_traj.gro  # The zero answers the first prompt given by trjconv. It
# corresponds to doing a conversion of the entire system. Delete echo and the pipe to see the other options.

Cylindricity_Traj.py -i wiggle_traj.gro