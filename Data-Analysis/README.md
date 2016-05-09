This folder contains scripts which analyze the results of simulations

The 'Cylindricity' scripts are named so because they were originally developed to calculate how close to a perfect cylinder the 
pore are after simulation. We had seen some clumping, micelle-like behavior for membrane structures having 5 monomer layers. The 
output of the script has deviated a bit from the original intention and includes the following analyses by tracking the postions of the sodium ions and carbonyl carbons:
  - Pore to Pore distance (based on central axis found by averaging x, y coordinates of atoms of interest - sodium or carbon)
  - Distribution of ions/carbonyl carbons with respect to the central axis
  - Density of ions/carbonyl carbons with respect to central axis

- Cylindricity.py is the original script and looks at the output .gro file
- Cylindricity_Traj.py looks at the entire simulation trajectory

Animated_Plot.py was made for Cylindricity_Traj.py and animates how the Pore-to-Pore distances change over the course of a 
trajectory. It is already incorporated into Cylindricty_Traj.py but it is a good script to understand if you want to make 
animated plots -- which can be a powerful tool

Single_monomer.py and Box_vectors.py are both scripts used to view outputs from GROMACS in the form of .xvg files. .xvg files are
output from gmx energy. 
