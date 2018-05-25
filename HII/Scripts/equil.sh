#!/bin/bash

set -e  # exit immediately after error

BUILD_MON="NAcarb11V"
start_config="initial.gro"
ring_restraints="C C1 C2 C3 C4 C5"
forces="3162 56 8 3 2 1 0"
MPI="off"
NP=4
T=300
equil_length=1000000  # equilibrium simulation length
python="python3"
quit_early=0
restraint_residue='HII'

while getopts "b:r:m:t:p:f:e:S:P:q:R:" opt; do
    case $opt in
    b) BUILD_MON=$OPTARG;;
    r) ring_restraints=$OPTARG;;
    m) MPI=$OPTARG;;
    t) T=$OPTARG;;
    p) NP=$OPTARG;;
    f) forces=$OPTARG;;
    e) equil_length=$OPTARG;;
    S) start_config=$OPTARG;;
    P) python=$OPTARG;;
    q) quit_early=$OPTARG;;
    R) restraint_residue=$OPTARG;;
    esac
done

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [[ ${MPI} == "on" ]]; then
    GMX="mpirun -np ${NP} gmx_mpi"
else
    GMX="gmx"
fi

# make input files
${python} ${DIR}/input.py -b ${BUILD_MON} -l 50 --restraints ${restraint_residue} --temp ${T} -f 50 --genvel yes -c ${start_config}
# create restrained topology

${python} ${DIR}/restrain.py -f 1000000 1000000 1000000 -A xyz -g ${start_config} -m ${BUILD_MON} -a ${ring_restraints}

gmx grompp -f em.mdp -p topol.top -c ${start_config} -o em
${GMX} mdrun -v -deffnm em
gmx grompp -f npt.mdp -p topol.top -c em.gro -o npt
${GMX} mdrun -v -deffnm npt

cp npt.gro 1000000.gro
cp npt.trr 1000000.trr

if [[ ${quit_early} -eq 1 ]]; then
    exit
fi

${python} ${DIR}/input.py -b ${BUILD_MON} -l 50 --restraints ${restraint_residue} --temp ${T} -f 50 --genvel no -c ${start_config}

for f in ${forces}; do

    ${python} ${DIR}/restrain.py -f ${f} ${f} ${f} -A xyz -g ${start_config} -m ${BUILD_MON} -a ${ring_restraints}
    gmx grompp -f npt.mdp -p topol.top -c npt.gro -o npt
    ${GMX} mdrun -v -deffnm npt

	cp npt.gro ${f}.gro
	cp npt.trr ${f}.trr

done

${python} ${DIR}/input.py -b ${BUILD_MON} -l 5000 --temp ${T} -f 50 --genvel no -c ${start_config} --barostat berendsen

gmx grompp -f npt.mdp -p topol.top -c npt.gro -o berendsen
${GMX} mdrun -v -deffnm berendsen
${python} ${DIR}/input.py -b ${BUILD_MON} -l ${equil_length} --temp ${T} -f 10000 --barostat Parrinello-Rahman --genvel no -c berendsen.gro
gmx grompp -f npt.mdp -p topol.top -c berendsen.gro -o PR
${GMX} mdrun -v -deffnm PR