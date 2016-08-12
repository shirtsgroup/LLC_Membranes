#!/bin/bash

#SBATCH -p GPU-shared --gres=gpu:1
#SBATCH -N 1 --ntasks-per-node 1
#SBATCH --time 00:30:00

module load gromacs
module load python/2.7.11_gcc

mpirun -np 1 gmx_mpi mdrun -v -deffnm 75box_emf_w1 -ntomp 7
