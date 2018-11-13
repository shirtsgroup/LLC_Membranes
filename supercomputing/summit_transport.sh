#!/bin/bash

# Add solutes to a solvated configuration then equilibrate the system
# If the system is cross-linked you will need to provide a cross-linked configuration (xlinked.gro) and an .itp topology file (assembly.itp)

#SBATCH -p sgpu
#SBATCH -N 1 
#SBATCH -t 24:00:00  # max wall time is 24 hrs

module load gromacs/2018_gpu
export GMX_MAXBACKUP=-1

place_solutes_pores.py -g xlinked.gro -o initial.gro -r ETH -n 6 -mdps -mpi 4  # replace -r ETH with appropriate solute

# fully energy minimize initial configuration
mpirun -np 1 gmx_mpi grompp -f em.mdp -p topol.top -c initial.gro -o initial_em
mpirun -np 4 gmx_mpi mdrun -v -deffnm initial_em

# run 5 ns berendsen simulation
mpirun -np 1 gmx_mpi grompp -f berendsen.mdp -p topol.top -c initial_em.gro -o berendsen
mpirun -np 4 gmx_mpi mdrun -v -deffnm berendsen

# run 400 ns Parinello-Rahman simulation
mpirun -np 1 gmx_mpi grompp -f PR.mdp -p topol.top -c berendsen.gro -o PR
mpirun -np 4 gmx_mpi mdrun -v -deffnm PR
