#!/bin/bash

#SBATCH -p GPU --gres=gpu:k80:4
#SBATCH -N 1 --tasks-per-node=4
#SBATCH -t 48:00:00

module load gromacs
module load gromacs/2018_gpu
module load python/3.6.4_gcc5_np1.14.5
