#!/bin/bash

set -e  # exit immediately after error

# default values
BUILD_MON="NAcarb11V"
start_config="initial.gro"
ring_restraints="C C1 C2 C3 C4 C5"  # change these
forces="3162 56 8 3 2 1 0"
MPI="off"
NP=4
T=300
equil_length=1000000  # equilibrium simulation length
python="python"
quit_early=0
restraint_residue='HII'
pd=0
ncol=5

OPTINID=1
while getopts "b:r:m:t:p:f:e:S:P:q:R:n:o:" opt; do
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
    n) ncol=$OPTARG;;  # number of columns per pore
    o) pd=$OPTARG;;
    esac
done

export GMX_MAXBACKUP=-1  # don't save backup files

# directory where this script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# If you want to use MPI, then the syntax for GROMACS changes slightly
if [[ ${MPI} == "on" ]]; then
    GMX="mpirun -np ${NP} gmx_mpi"
else
    GMX="gmx"
fi

# assumes defaults in build.py are okay. Adjust them if needed
#${python} ${DIR}/build.py -pd ${pd} -b ${BUILD_MON}.gro -c ${ncol}

# make input files
${python} ${DIR}/input.py -b ${BUILD_MON} -l 50 --restraints ${restraint_residue} --temp ${T} -f 50 --genvel yes -c ${start_config} -s 50000

# create restrained topology
${python} ${DIR}/restrain.py -f 1000000 1000000 1000000 -A xyz -g ${start_config} -m ${BUILD_MON} -a ${ring_restraints}

# try energy minimization
gmx grompp -f em.mdp -p topol.top -c ${start_config} -o em -r ${start_config}
${GMX} mdrun -v -deffnm em

# check for negative potential energy after minimization. If it's positive, restart the script.
n=$(awk '/Potential Energy/ {print $4}' em.log)
echo $n
if [[ ${n:0:2} == *"."* ]]; then
    #bash ${DIR}/equil.sh -b ${BUILD_MON} -r ${ring_restraints} -m ${MPI} -t ${T} -p ${NP} -f ${forces} -e ${equil_length} -S ${start_config} -P ${python} -q ${quit_early} -R ${restraint_residue} -n ${ncol} -o ${pd}
    bash ${DIR}/equil.sh -m ${MPI} -p ${NP} -P ${python} -q ${quit_early} -n ${ncol} -o ${pd}
fi

# run 1st restrained equilibration
gmx grompp -f npt.mdp -p topol.top -c em.gro -o npt -r em.gro
${GMX} mdrun -v -deffnm npt

# change name to reflect position restraint force constant
cp npt.gro 1000000.gro
cp npt.trr 1000000.trr

# if you want, stop the equilibration just after this step
if [[ ${quit_early} -eq 1 ]]; then
    exit
fi

# change .mdp file so it does not generate velocities
${python} ${DIR}/input.py -b ${BUILD_MON} -l 50 --restraints ${restraint_residue} --temp ${T} -f 50 --genvel no -c ${start_config}

# gradually reduce force constants after 50 ps increments
for f in ${forces}; do

    ${python} ${DIR}/restrain.py -f ${f} ${f} ${f} -A xyz -g ${start_config} -m ${BUILD_MON} -a ${ring_restraints}
    gmx grompp -f npt.mdp -p topol.top -c npt.gro -o npt -r npt.gro
    ${GMX} mdrun -v -deffnm npt

	cp npt.gro ${f}.gro
	cp npt.trr ${f}.trr

done

# run berendsen equilibration
${python} ${DIR}/input.py -b ${BUILD_MON} -l 5000 --temp ${T} -f 50 --genvel no -c ${start_config} --barostat berendsen
gmx grompp -f npt.mdp -p topol.top -c npt.gro -o berendsen
${GMX} mdrun -v -deffnm berendsen

# run Parinello-Rahman equilibration
${python} ${DIR}/input.py -b ${BUILD_MON} -l ${equil_length} --temp ${T} -f 10000 --barostat Parrinello-Rahman --genvel no -c berendsen.gro
gmx grompp -f npt.mdp -p topol.top -c berendsen.gro -o PR
${GMX} mdrun -v -deffnm PR
