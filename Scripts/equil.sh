#!/bin/bash

set -e  # exit immediately after error

#restrain.py -f 1000000 -A xyz -r on -D off -w off -g initial.gro --novsites -m NAcarb11V
#gmx grompp -f em.mdp -p topol.top -c box.gro -o em
#gmx mdrun -v -deffnm em
restrain.py -f 1000000 -A xyz -r on -D off -w off --novsite -m NAcarb11V
gmx grompp -f npt.mdp -p topol.top -c em.gro -o npt
gmx mdrun -v -deffnm npt
cp npt.gro 1000000.gro
cp npt.trr 1000000.trr

for f in 3162 56 8 3 2 1 0; do
	restrain.py -f ${f} -A xyz -r on -D off -w off --novsites -m NAcarb11V
	gmx grompp -f npt.mdp -p topol.top -c npt.gro -o npt
	gmx mdrun -v -deffnm npt
	cp npt.gro ${f}.gro
	cp npt.trr ${f}.trr
done
