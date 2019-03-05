#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # the directory where this script is located.
DIR2="${DIR}/../../HII/Scripts"  # some scripts aren't located in the same directory as $DIR

build_mon="Dibrpyr14"  # monomer to build system with
phase="gyroid"  # bcc phase to build
dimension=10  # length of cubic unitcell box vector
density=1.1  # density g/cm^3
SHIFT=0  # shift head group position along normal (lowercase 'shift' is bash built-in command)
ion="BR"
wt_percent=77.1  # weight percent of monomer
solvent='glycerol'  # name of solvent with bcc monomer
top='topol.top'  # name of topology
temp=343
restraint_atoms="C23 C17 C28"
sol_res=GLY # solvent residue name
cluster=0  # 0 if run without mpi, 1 if run with mpi (regular MPI or with GPU)
NP=4  # number of MPI processes
python='python3'

while getopts "b:p:d:r:n:I:w:S:T:R:c:P:" opt; do
    case $opt in
    b) build_mon=$OPTARG;;
    p) phase=$OPTARG;;
    d) dimension=$OPTARG;;
    r) density=$OPTARG;;
    n) SHIFT=$OPTARG;;
    I) ion=$OPTARG;;
    w) wt_percent=$OPTARG;;
    S) solvent=$OPTARG;;
    T) temp=$OPTARG;;
    R) restraint_atoms=$OPTARG;;
    c) cluster=$OPTARG;;
    P) NP=$OPTARG;;
    esac
done

echo $NP
echo ${cluster}
if [[ ${cluster} -eq 1 ]]; then
    GMX="mpirun -np ${NP} gmx_mpi"
else
    GMX="gmx"
fi

echo $GMX
echo $NP
exec ./test.sh -c 1 -P 4

exit
{GMX} --version 
