## How to parameterize an arbitrary monomer with GAFF using antechamber ##

1. Create structure in ChemDraw (or any other tool i.e. Avogadro). Save as *.mol (choose MDL molfile, not MDL molfile V3000)
2. Convert to .pdb using <a href="https://cactus.nci.nih.gov/translate/">online smiles converter</a>
3. Make sure the .pdb has the residue name in the correct location (starts at character 18 in ATOM section). 
    1. If it is blank, you will get a difficult-to-trace memory allocation error like "_memory allocation error for *bond_array_"
4. Parameterize with antechamber following the <a href="http://ambermd.org/tutorials/basic/tutorial4b/">tutorial</a>.
5. Convert the output *.prmtop and *.inpcrd files to gromacs *.gro and *.top file using acpype.py
6. Energy minimze the new structure
_NOTE: Steps 4 - 6 can be done automatically using param.sh_
