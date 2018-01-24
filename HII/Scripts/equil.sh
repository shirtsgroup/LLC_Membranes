#!/bin/bash

set -e  # exit immediately after error

BUILD_MON="NAcarb11V"
start_config="initial.gro"
x=9.5
y=9.5
z=7.4
ring_restraints="C C1 C2 C3 C4 C5"
forces="3162 56 8 3 2 1 0"
MPI="off"
NP=4
T=300
equil_length=1000000  # equilibrium simulation length
solvate=0
python="python3"

while getopts "b:x:y:z:r:m:t:p:f:e:s:S:" opt; do
    case $opt in
    b) BUILD_MON=$OPTARG;;
    x) x=$OPTARG;;
    y) y=$OPTARG;;
    z) z=$OPTARG;;
    r) ring_restraints=$OPTARG;;
    m) MPI=$OPTARG;;
    t) T=$OPTARG;;
    p) NP=$OPTARG;;
    f) forces=$OPTARG;;
    e) equil_length=$OPTARG;;
    s) solvate=$OPTARG;;
    S) start_config=$OPTARG;;
    esac
done

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [[ $solvate -eq 1 ]]; then
    ${python} ${DIR}/input.py -b ${BUILD_MON} -l 50 --restraints --temp ${T} -f 50 --genvel yes -S -c ${start_config}
else
    ${python} ${DIR}/input.py -b ${BUILD_MON} -l 50 --restraints --temp ${T} -f 50 --genvel yes -c ${start_config}
fi

${python} ${DIR}/restrain.py -f 1000000 -A xyz -r on -D off -w off -g ${start_config} --novsites -m ${BUILD_MON} -a ${ring_restraints}

if [ "${MPI}" == "on" ]; then
    gmx editconf -f ${start_config} -o box.gro -c -bt triclinic -box ${x} ${y} ${z} -angles 90 90 120
    gmx grompp -f em.mdp -p topol.top -c box.gro -o em
    mpirun -np ${NP} gmx_mpi mdrun -v -deffnm em
    gmx grompp -f npt.mdp -p topol.top -c em.gro -o npt
    mpirun -np ${NP} gmx_mpi mdrun -v -deffnm npt
else
    gmx editconf -f ${start_config} -o box.gro -c -bt triclinic -box ${x} ${y} ${z} -angles 90 90 120
    gmx grompp -f em.mdp -p topol.top -c box.gro -o em
    gmx mdrun -v -deffnm em
    gmx grompp -f npt.mdp -p topol.top -c em.gro -o npt
    gmx mdrun -v -deffnm npt
fi

cp npt.gro 1000000.gro
cp npt.trr 1000000.trr

${python} ${DIR}/input.py --mdp -b ${BUILD_MON} -l 50 --restraints --temp ${T} -f 50 --genvel no  -c ${start_config} # use velocities from previous sim

for f in ${forces}; do
	${python} ${DIR}/restrain.py -f ${f} -A xyz -r on -D off -w off --novsites -m ${BUILD_MON} -a ${ring_restraints} -g ${start_config}
	if [ ${MPI} == 'on' ]; then
            gmx grompp -f npt.mdp -p topol.top -c npt.gro -o npt
	    mpirun -np ${NP} gmx_mpi mdrun -v -deffnm npt
	else
	    gmx grompp -f npt.mdp -p topol.top -c npt.gro -o npt
	    gmx mdrun -v -deffnm npt
	fi
	cp npt.gro ${f}.gro
	cp npt.trr ${f}.trr
done

#input.py --mdp -b ${BUILD_MON} -l 20000 --restraints --temp ${T} -f 50 --genvel no  # use velocities from previous sim
#gmx grompp -f npt.mdp -p topol.top -c npt.gro -o wiggle
if [[ $solvate -eq 1 ]]; then
    ${python} ${DIR}/input.py -b ${BUILD_MON} -l 5000 --temp ${T} -f 50 --barostat berendsen --genvel no -S -c ${start_config}
else
    ${python} ${DIR}/input.py -b ${BUILD_MON} -l 5000 --temp ${T} -f 50 --barostat berendsen --genvel no -c ${start_config}
fi

if [ ${MPI} == 'on' ]; then
    gmx grompp -f npt.mdp -p topol.top -c npt.gro -o berendsen
    mpirun -np ${NP} gmx_mpi mdrun -v -deffnm berendsen
    ${python} ${DIR}/input.py -b ${BUILD_MON} -l ${equil_length} --temp ${T} -f 10000 --barostat Parrinello-Rahman --genvel no
    gmx grompp -f npt.mdp -p topol.top -c berendsen.gro -o PR
else
    gmx grompp -f npt.mdp -p topol.top -c 0.gro -o berendsen # run it out for a bit
    gmx mdrun -v -deffnm berendsen
    ${python} ${DIR}/input.py -b ${BUILD_MON} -l ${equil_length} --temp ${T} -f 10000 --barostat Parrinello-Rahman --genvel no
    gmx grompp -f npt.mdp -p topol.top -c berendsen.gro -o PR
fi

#input.py -b ${BUILD_MON} -l ${equil_length} --temp ${T} -f 1000 --barostat Parrinello-Rahman --genvel no
#gmx grompp -f npt.mdp -p topol.top -c wiggle.gro -o wiggle

#analysis.sh # get pore spacing, thickness, pore size, convert the trajectory to be used in XrayDiffraction.exe

