#!/bin/bash

# Use the bits of this that are required for whatever you are doing. The following is capable of running python scripts from this repository as well as gromacs compiled for threaded and GPU runs.

#SBATCH -p GPU --gres=gpu:k80:4
#SBATCH -N 1 --tasks-per-node=4
#SBATCH -t 48:00:00

module purge  # clear out modules or else mdtraj gets screwed up
module load gromacs/2018_gpu  # GPU accelerated gromacs installed by Marcela Madrid
module load python/3.6.4_gcc5_np1.14.5  # this is the only version where numpy works. Will need to install mdtraj (pip install --user mdtraj)
export PYTHONPATH="$PYTHONPATH:/home/bjc/Gromacs"  # Bridges forgets what is written in your .bashrc
source /home/bjc/pkgs/gromacs/2018.3/bin/GMXRC # User installed threaded version of Gromacs since running grompp in parallel probably isn't a good idea
