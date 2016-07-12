#! /bin/bash

echo "#SBATCH --qos janus-long" > slurm.sh
echo "#SBATCH --nodes 8" >> slurm.sh
echo "#SBATCH --ntasks-per-node 1" >> slurm.sh
echo "#SBATCH --time 168:00:00" >> slurm.sh

