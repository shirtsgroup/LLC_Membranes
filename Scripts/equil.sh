#!/bin/bash

set -e  # exit immediately after error

restrain.py -f 100000 -A xyz -r on -D off -w off -g initial.gro
gmx grompp -f em.mdp -p topol.top -c box.gro -o em
gmx mdrun -v -deffnm em
restrain.py -f 100000 -A xyz -r on -D off -w off
gmx grompp -f npt.mdp -p topol.top -c em.gro -o npt
gmx mdrun -v -deffnm npt
cp npt.gro 100000.gro
cp npt.trr 100000.trr

for f in 3162 56 8 3 2 1 0; do
	restrain.py -f ${f} -A xyz -r on -D off -w off
	gmx grompp -f npt.mdp -p topol.top -c npt.gro -o npt
	gmx mdrun -v -deffnm npt
	cp npt.gro ${f}.gro
	cp npt.trr ${f}.trr
done