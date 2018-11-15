#!/bin/bash

#SBATCH -p GPU --gres=gpu:k80:4
#SBATCH -N 1 --tasks-per-node=4
#SBATCH -t 48:00:00

module load gromacs/2018_gpu
mpirun -np 4 gmx_mpi mdrun -s PR.tpr -cpi PR.cpt -append -v -deffnm PR
