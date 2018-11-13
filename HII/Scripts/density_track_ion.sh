#!/usr/bin/env bash

while getopts "p:" opt; do
    case $opt in
    p) PATH=$OPTARG;;
    esac
done
scp bjc@bridges.psc.edu:$PATH/\{$FILE1,$FILE2,$FILE3\} ./
#scp bjc@bridges.psc.edu:/pylon1/ct4s8bp/bjc/Solvate/${PATH}/\{wiggle.tpr,wiggle.trr,wiggle.gro\} ./


~/PycharmProjects/GitHub/Scripts/Track_Ion.py -t wiggle.trr -s wiggle.tpr -g wiggle.gro

~/PycharmProjects/GitHub/Scripts/density.py -t wiggle.trr -s wiggle.tpr -g wiggle.gro -b buffer
