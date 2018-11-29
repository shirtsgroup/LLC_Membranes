#!/bin/bash

# Equilibrate a hexagonal columnar phase liquid crystal without water

#SBATCH -p sgpu  # use the GPU partition
#SBATCH -N 1  # number of nodes to use
#SBATCH -t 24:00:00  # max wall time is 24 hrs

module load gromacs/2018_gpu

equil.py -b monomer.gro -mpi -np 4 --ring_restraints C C1 C2 C3 C4 C5 --restraint_residue HII
