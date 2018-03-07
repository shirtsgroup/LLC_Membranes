#!/usr/bin/env bash

set -e # exit on error

nmon=5  # number of monomers per layer
radius=8  # radius of pore (angstroms)
offset=1  # 1 for parallel displaced, else for sandwiched
gap_size=3  # size of water gap (nm)
nwater_tails=600  # number of waters to add to tail region
pore_water=1.09  # weight percent of water to pores
gap=0  # set to 1 if you want a water gap

while getopts "r:m:o:g:t:p:G:" opt; do
    case $opt in
    r) radius=$OPTARG;;
    m) nmon=$OPTARG;;
    o) offset=$OPTARG;;
    g) gap_size=$OPTARG;;
    t) nwater_tails=$OPTARG;;
    p) pore_water=$OPTARG;;
    G) gap=$OPTARG;;
    esac
done

####################################  SOLVATION PROCEDURE  ############################################

# create initial structure
if [[ ${offset} -eq 1 ]]; then
    build.py -m ${nmon} -r ${radius} -O
else
    build.py -m ${nmon} -r ${radius}
fi

# run short restrained equilibration to loosen up tails
equil.sh -q 1 -x 9 -y 9 -z 7.4

if [[ ${gap} -eq 1 ]]; then
    # add gap to initial structure but don't solvate
    add_gap.py -g 1000000.gro -l ${gap_size}
    # solvate pores with gap. This will create a water bath and a selectively controlled initial water concentration within the pores
    solvate_pore.py -g gap.gro -l ${gap_size} -o solvated_pore.gro -r ${radius} -w ${pore_water}
else
    solvate_pore.py -g 1000000.gro -o solvated_pore.gro -r ${radius} -w ${pore_water}
fi

# make input files again so that water is included in the topology appropriately
input.py -c solvated_pore.gro -S

gmx grompp -f em.mdp -p topol.top -c solvated_pore.gro -o em_solvated_pore
gmx mdrun -v -deffnm em_solvated_pore

# solvate tails
solvate_tails.py -g em_solvated_pore.gro -p topol.top -o fully_solvated.gro -n ${nwater_tails}

# recreated topology and set up npt.mdp for a 5 ns berendsen pressure controlled simulation
input.py -c fully_solvated.gro -S -l 5000

# energy minimize
gmx grompp -f em.mdp -p topol.top -c fully_solvated.gro -o em_fully_solvated
gmx mdrun -v -deffnm em_fully_solvated

# Set up 5 ns berendsen controlled run (these parameters were set by default earlier) and 400 ns Parrinello-Rahman run
gmx grompp -f npt.mdp -p topol.top -c em_fully_solvated.gro -o berendsen
input.py -c fully_solvated.gro -S -l 400000 --barostat Parrinello-Rahman --genvel no -f 10000

exit
# Transfer to bridges to run it
scp berendsen.tpr npt.mdp topol.top bjc@bridges.psc.edu:/pylon5/ct4s8bp/bjc/Solvate/pores_tails/${radius}
ssh bjc@bridges.psc.edu "sbatch /pylon5/ct4s8bp/bjc/Solvate/pores_tails/${radius}/berendsen.sh"