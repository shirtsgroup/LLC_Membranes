This folder contain Shell scripts mostly made for the purpose of running simulations as simply and customizable as possible.

* Please note that you may need to modify some of the paths in the folders to match your computer's set-up

Anything that does not end in .sh are files that will be used in the shell scripts:

- wiggle.mdp      : .mdp file for simulations of pore assembly in vacuum
- em.mdp          : .mdp file for energy minimzation
- wiggle_solv.mdp : .mdp file for simulations of pore assembly solvated in water 
- NaPore.top      : topology file for vacuum pore assembly
- NaPore_water.top: topology file for solvated system (has water force field enabled)
- vdwradii.dat    : modified van der waals radius file so that water is not placed in places it shouldn't go when solvating

Quick guide to shell script name abbrevations:

- BS        : Build and Simulate
- BSS       : Build, Solvate and Simulate
- PJV       : Prepare Janus Input File in Vacuum (A truncation of BS)
- PJS       : Prepare Janus Input File for Solvated System (A truncation of BSS)
- ifv.sh    : Input files needed for a vacuum system
- ifs.sh    : Input files needed for a solvated system
- janus.sh  : ssh into janus by typing this instead of having to write the address every time
- bridges.sh: ssh into bridges

Description of Shell Scripts:

BS: Use to simulate a vacuum system. Builds structure from specified monomer, energy minimizes it, then runs a simulation for a duration specified by wiggle.mdp

BSS: Use to simulate a solvated system. Build structure from specified monomer, energy minimizes, runs short simulation (as defined by wiggle.mdp), solvates the system, energy minimizes again, runs a simulation with the solvated system.

PJV: Same as BS except the simulation is not run at the end. Instead, an .tpr file is output which can be fed to gromacs later

PJS: Same as BSS except the solvated simulation is not run. A .tpr file is left so that the simulation can be run on Janus

ifv.sh: Puts all of the files necessary for a vacuum simulation set-up into the current working directory

ifs.sh: Puts all of the files necessary for a solvated system simulation set-up into the current working directory

janus.sh/bridges.sh : one-liners that make is so you don't have to type ssh username@rc.colorado.edu -- it gets annoying
