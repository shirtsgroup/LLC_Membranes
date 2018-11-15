#!/bin/bash

# Add solutes to a solvated configuration then equilibrate the system
# If the system is cross-linked you will need to provide a cross-linked configuration (xlinked.gro) and an .itp topology file (assembly.itp)

#SBATCH -p GPU --gres=gpu:k80:4
#SBATCH -N 1 --tasks-per-node=4
#SBATCH -t 48:00:00

module load gromacs/2018_gpu  # GPU accelerated gromacs installed by Marcela Madrid
module load python/3.6.4_gcc5_np1.14.5  # this is the only version where numpy works. Will need to install mdtraj (pip install --user mdtraj)
export PYTHONPATH="$PYTHONPATH:/home/bjc/Gromacs"  # Bridges forgets what is written in your .bashrc
export GMX_MAXBACKUP=-1  # do not make back-ups

place_solutes_pores.py -g xlinked.gro -o initial.gro -r ETH -n 6 -mdps -mpi 4  # will parallelize poorly but MPI mdoule messes up python

module load mpi/intel_mpi  # loaded after python 
mpirun -np 1 gmx_mpi grompp -f berendsen.mdp -p topol.top -c initial.gro -o berendsen
mpirun -np 4 gmx_mpi mdrun -v -deffnm berendsen

mpirun -np 1 gmx_mpi grompp -f PR.mdp -p topol.top -c berendsen.gro -o PR
mpirun -np 4 gmx_mpi mdrun -v -deffnm P
