#!/bin/bash
set -e # exit upon error

export AMBERHOME=/home/NormaKLangdon/Programs/amber16
export PATH="$HOME/Programs/amber16/bin/:$PATH"
export PATH="$HOME/Programs/amber16/AmberTools/src/leap/:$PATH"

antechamber -i monomer.pdb -fi pdb -o monomer.mol2 -fo mol2 -c bcc -s 2 -nc -1 
parmchk -i monomer.mol2 -f mol2 -o monomer.frcmod

tleap -f tleap.in

./acpype.py -p monomer.prmtop -x monomer.inpcrd

python topedit.py

gmx editconf -f MOL_GMX.gro -o monB.gro -c -d 3 -bt cubic
gmx grompp -f Emin.mdp -p MOL_GMX.top -c monB.gro -o monBE.tpr
gmx mdrun -v -deffnm monBE
gmx grompp -f anneal.mdp -p MOL_GMX.top -c monBE.gro -o monBEA.tpr
gmx mdrun -v -deffnm monBEA

#Choose trajectory you want
