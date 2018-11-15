#!/bin/bash

#SBATCH -p sgpu
#SBATCH -N 1
#SBATCH -t 24:00:00

module load gromacs/2018_gpu
mpirun -np 4 gmx_mpi mdrun -s PR.tpr -cpi PR.cpt -append -v -deffnm PR
