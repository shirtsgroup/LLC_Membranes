#! /bin/bash

# Parameterize and energy minimize a monomer with GAFF using antechamber
# given a .pdb file 

set -e # exit upon error

name='monomer' # name of molecule. This name will carry through to output
nc=0  # net charge
res='MON' # name of residue being made, NOTE: This must match the residue in the .pdb file or you will get errors
anneal='no'  # change to 'yes' if you want a thermal annealing process carried out after energy minimization
input_path='.'

while getopts "n:c:r:a:p:" opt; do
    case $opt in
    n) name=$OPTARG;;
    c) nc=$OPTARG;;
    r) res=$OPTARG;;
    a) anneal=$OPTARG;;
    p) input_path=$OPTARG;;
    esac
done

antechamber -i ${name}.pdb -fi pdb -o ${name}.mol2 -fo mol2 -c bcc -s 2 -nc ${nc}  # The .pdb must have connectivity info!
# -c bcc tells antechamber to use AM1-BCC charge model
# -s flag just defines verbosity of output
parmchk -i ${name}.mol2 -f mol2 -o ${name}.frcmod

# Create input to tleap
echo "source oldff/leaprc.ff99SB" > tleap.in
echo "source leaprc.gaff" >> tleap.in  # make sure tleap knows about GAFF forcefield
echo "${res} = loadmol2 ${name}.mol2" >> tleap.in  # load monomer
echo "check ${res}" >> tleap.in  # checks for missing parameters
echo "loadamberparams ${name}.frcmod" >> tleap.in  # tell tLeap what to do with missing parameters
echo "saveoff ${res} ${res}.lib" >> tleap.in  # create a library file for the new residue
echo "saveamberparm ${res} ${name}.prmtop ${name}.inpcrd" >> tleap.in  # make topology files
echo "quit" >> tleap.in  # exit out of tleap

tleap -f tleap.in  # run the previous block in tleap

#sed -i -e "s/${res: -3}/${res}/g" ${name}.prmtop
#sed -i -e "s/${res: -3}/${res}/g" ${name}.inpcrd

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # the directory where this script is located.
DIR=${input_path}

acpype.py -p ${name}.prmtop -x ${name}.inpcrd

# Rename. I prefer to get rid of the GMX part. Copy them though so we have backup
cp ${res}_GMX.gro ${res}.gro
cp ${res}_GMX.top ${res}.top

gmx editconf -f ${res}.gro -o box.gro -c -d 3 -bt cubic  # create a cubic box with lots of space for the monomer
gmx grompp -f ${DIR}/em.mdp -p ${res}.top -c box.gro -o em  # create atomic level input file
gmx mdrun -v -deffnm em  # run energy minimization
gmx editconf -f em.gro -o ${res}_new.pdb  # will need this for molcharge later

if [ "${anneal}" == "on" ]; then  # anneal if necessary
	gmx grompp -f ${DIR}/anneal.mdp -p ${res}.top -c em.gro -o anneal.tpr
	gmx mdrun -v -deffnm anneal
	gmx editconf -f anneal.gro -o ${res}_new.pdb  # will need this for molcharge later
fi

cp ${name}.mol2 ${name}_bak.mol2  # backup mol2 initially generated in case the names are the same

# Use OpenEye to assign new partial charges
molcharge -in ${res}_new.pdb -out ${res}_new.mol2 -method am1bccsym

# replace charges output by antechamber with those output by molcharge
insertmol2charges.py -m ${res}_new.mol2 -t ${res}.top -c ${nc} -o ${res}.top

# Energy minimize again
gmx grompp -f ${DIR}/em.mdp -p ${res}.top -c ${res}_new.pdb -o em
gmx mdrun -v -deffnm em

cp em.gro ${res}.gro
cp ${res}.gro ${name}.gro
cp ${res}.top ${name}.top