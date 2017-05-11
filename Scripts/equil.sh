#!/bin/bash

set -e  # exit immediately after error

BUILD_MON="NAcarb11V"
x=9
y=9
z=8
ring_restraints="C C1 C2 C3 C4 C5"
MPI="off"
NP=4
T=300

while getopts "b:x:y:z:r:m:t:p:" opt; do
    case $opt in
    b) BUILD_MON=$OPTARG;;
    x) x=$OPTARG;;
    y) y=$OPTARG;;
    z) z=$OPTARG;;
    r) ring_restraints=$OPTARG;;
    m) MPI=$OPTARG;;
    t) T=$OPTARG;;
    p) NP=$OPTARG;;
    esac
done

input.py -b ${BUILD_MON} -l 50 --restraints --temp ${T} -f 50 --genvel yes
restrain.py -f 1000000 -A xyz -r on -D off -w off -g initial.gro --novsites -m ${BUILD_MON} -a ${ring_restraints}

if [ "${MPI}" == "on" ]; then
    gmx_mpi editconf -f initial.gro -o box.gro -c -bt triclinic -box ${x} ${y} ${z} -angles 90 90 120
    gmx_mpi grompp -f em.mdp -p topol.top -c box.gro -o em
    mpirun -np ${NP} gmx_mpi mdrun -v -deffnm em
    gmx_mpi grompp -f npt.mdp -p topol.top -c em.gro -o npt
    mpirun -np ${NP} gmx_mpi mdrun -v -deffnm npt
else
    gmx editconf -f initial.gro -o box.gro -c -bt triclinic -box ${x} ${y} ${z} -angles 90 90 120
    gmx grompp -f em.mdp -p topol.top -c box.gro -o em
    gmx mdrun -v -deffnm em
    gmx grompp -f npt.mdp -p topol.top -c em.gro -o npt
    gmx mdrun -v -deffnm npt
fi

cp npt.gro 1000000.gro
cp npt.trr 1000000.trr

input.py --mdp -b ${BUILD_MON} -l 50 --restraints --temp ${T} -f 50 --genvel no  # use velocities from previous sim

for f in 3162 56 8 3 2 1 0; do
	restrain.py -f ${f} -A xyz -r on -D off -w off --novsites -m ${BUILD_MON} -a ${ring_restraints}
	if [ ${MPI} == 'on' ]; then
		gmx_mpi grompp -f npt.mdp -p topol.top -c npt.gro -o npt
	    mpirun -np ${NP} gmx_mpi mdrun -v -deffnm npt
	else
	    gmx grompp -f npt.mdp -p topol.top -c npt.gro -o npt
	    gmx mdrun -v -deffnm npt
	fi
	cp npt.gro ${f}.gro
	cp npt.trr ${f}.trr
done

input.py -b ${BUILD_MON} -l 50 --temp ${T} -f 50 --barostat berendsen --genvel no # put pressure control back on
gmx grompp -f npt.mdp -p topol.top -c 0.gro -o wiggle # run it out for a bit
gmx mdrun -v -deffnm wiggle

input.py -b ${BUILD_MON} -l 100000 --temp ${T} -f 100 --barostat Parrinello-Rahman --genvel no
gmx grompp -f npt.mdp -p topol.top -c wiggle.gro -o wiggle
#analysis.sh # get pore spacing, thickness, pore size, convert the trajectory to be used in XrayDiffraction.exe

