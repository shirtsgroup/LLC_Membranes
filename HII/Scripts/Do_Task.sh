#!/bin/bash

# When I want to run lots of things overnight but don't want to keep checking on them, write a script here:

# Run antechamber and parameterize 5 different monomers with GAFF:

cd /home/bcoscia/Documents/Parameterization/Monomer_configs

for dir in 1 2 3 4 5 6; do
    cd $dir
    antechamber -i monomer${dir}.pdb -fi pdb -o monomer${dir}.mol2 -fo mol2 -c bcc -s 2 -nc -1
    parmchk -i monomer${dir}.mol2 -f mol2 -o monomer${dir}.frcmod
    tleap -f tleap.in
    acpype.py -p monomer${dir}.prmtop -x monomer${dir}.inpcrd
    cd ..
done

cd /home/bcoscia/PycharmProjects/Github

mkdir 5_0_10

cd 5_0_10

cp ../Crosslink/wiggle_init.gro .
mv wiggle_init.gro wiggle.gro

xlink.sh -c 5 -t 0 -d 15

cd ../

mkdir 5_10_10

cd 5_10_10

cp ../Crosslink/wiggle_init.gro .
mv wiggle_init.gro wiggle.gro

xlink.sh -c 5 -t 10 -d 10

cd ../

mkdir 10_5_15

cd 10_5_15

cp ../Crosslink/wiggle_init.gro .
mv wiggle_init.gro wiggle.gro

xlink.sh -c 10 -t 5 -d 15