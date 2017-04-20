## How to parameterize an arbitrary monomer with GAFF using antechamber ##

1. Create structure in ChemDraw (or any other tool i.e. Avogadro). Save as *.mol (choose MDL molfile, not MDL molfile V3000)
2. Convert to .pdb using <a href="https://cactus.nci.nih.gov/translate/">online smiles converter</a>
3. *Make sure the .pdb has the residue name in the correct location (starts at character 18). 
       *If it is blank, you will get a difficult-to-trace memory allocation error like memory allocation error for *bond_array
4. Parameterize with antechamber following the <a href="http://ambermd.org/tutorials/basic/tutorial4b/">tutorial</a>.
