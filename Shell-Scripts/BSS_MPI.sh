#!/bin/bash
# Script to build a structure to a specified number of layers using a specified
# monomer, then energy minimize the structure, and finally run a simulation
# based on inputs to the .mdp files

# Copy in necessary files
ifs.sh

# Define Variables with Default values

# Choose which monomer to build with
MONOMER='monomer4.pdb'  # Structure file to be used

#MPI Options
NODES=16

# Energy minimization parameters:
INTEGRATOR_EM=steep  # Integrator for energy minimization
NSTEPS_EM=50000  # Maximum number of steps to take for energy minimization
CUTOFF_EM='verlet'  # Cut-off Scheme
NSTLIST=10  # Neighborlist - changed automatically by gromacs unless it is set equal to 1

# Membrane Dimensions
NO_MONOMERS=6  # Number of monomers in 1 layer
RADIUS=3  # Initial pore radius, angstroms
PORE2PORE=40  # Pore-to-Pore distance, angstroms
NOPORES=4  # Number of pores to be built
DBWL=10  # Distance Between Layers
LAYERS=20    # Number of layers wanted in the structure

# Box Vector Parameters
XVECT=8.0  # Box vector in the x direction
YVECT=8.0  # Box vector in the y direction (not this will be multiplied by sin(120) for a monoclinic cell of 120 degrees
INCREMENT=0.1  # Increment to increase the box vector by if there is a LINCS error

# Simulation Parameters
SIM_TITLE='Equilibration in Water'  # Title of simulation
CUTOFF_MD='verlet'  # Cut-off scheme for simulation
INTEGRATOR_MD='md'  # Integrator type for simulation
DT=0.002  # Time step (ps)
SIM_LENGTH=0.01  #length of intermediate 'wiggle' simulation (ns)
FRAMES=50  # Number of frames in trajectory
TCOUPL='v-rescale'
REF_T=300  # Reference Temperature, K
PCOUPL='berendsen'
PTYPE='semiisotropic'
REF_P=1  # Reference Pressure, bar
COMPRESSIBILITY=4.5e-5  # Isothermal compressibility, bar^-1
PBC='xyz'

# Solvation
WATER_LAYER=6  # thickness (nm) between membrane layers in the z direction
SOLV_LENGTH=1  # length of solvated simulation, nanoseconds

# Reference for flags associated with each variable

# -n  :   NODES ... number of nodes to use
# -M  :   MONOMER ... Structure file to be used
# -I  :   INTEGRATOR_EM ... Integrator for energy minimization
# -s  :   NSTEPS_EM ... Maximum number of steps to take for energy minimization
# -c  :   CUTOFF_EM ... Cut-off Scheme
# -t  :   NSTLIST ... Neighborlist - changed automatically by gromacs unless it is set equal to 1
# -o  :   NO_MONOMERS ... Number of monomers in 1 layer
# -r  :   RADIUS ... Initial pore radius, angstroms
# -p  :   PORE2PORE ... Pore-to-Pore distance, angstroms
# -P  :   NOPORES ... Number of pores to be built
# -w  :   DBWL ... Distance Between Layers
# -l  :   LAYERS ... Number of layers wanted in the structure
# -x  :   XVECT ... Box vector in the x direction
# -y  :   YVECT ... Box vector in the y direction (not this will be multiplied by sin(120) for a monoclinic cell of 120 degrees
# -e  :   INCREMENT ... Increment to increase the box vector by if there is a LINCS error
# -T  :   SIM_TITLE ... Title of simulation
# -C  :   CUTOFF_MD ... Cut-off scheme for simulation
# -i  :   INTEGRATOR_MD ... Integrator type for simulation
# -D  :   DT ... Time step (ps)
# -L  :   SIM_LENGTH ... Simulation length, nanoseconds (Vacuum)
# -f  :   FRAMES ... Number of frames in trajectory
# -v  :   TCOUPL ... Temperature Coupling
# -K  :   REF_T ... Reference Temperature, K
# -b  :   PCOUPL ... Pressure Coupling
# -Y  :   PTYPE ... i.e. Isotropic, semiisotropic
# -B  :   REF_P ... Reference Pressure, bar
# -R  :   COMPRESSIBILITY ... Isothermal compressibility, bar^-1
# -Z  :   PBC ... Periodic Boundary directions
# -V  :   SOLV_LENGTH ... length of final simulation

while getopts "M:I:S:c:t:o:r:p:P:w:l:x:y:e:T:C:i:D:L:f:v:K:b:Y:B:R:Z:V:" opt; do
    case $opt in
    n)  NODES=$OPTARG;;
    M)  MONOMER=$OPTARG;;
    I)  INTEGRATOR_EM=$OPTARG;;
    S)  NSTEPS_EM=$OPTARG;;
    c)  CUTOFF_EM=$OPTARG;;
    t)  NSTLIST=$OPTARG;;
    o)  NO_MONOMERS=$OPTARG;;
    r)  RADIUS=$OPTARG;;
    p)  PORE2PORE=$OPTARG;;
    P)  NOPORES=$OPTARG;;
    w)  DBWL=$OPTARG;;
    l)  LAYERS=$OPTARG;;
    x)  XVECT=$OPTARG;;
    y)  YVECT=$OPTARG;;
    e)  INCREMENT=$OPTARG;;
    T)  SIM_TITLE=$OPTARG;;
    C)  CUTOFF_MD=$OPTARG;;
    i)  INTEGRATOR_MD=$OPTARG;;
    D)  DT=$OPTARG;;
    L)  SIM_LENGTH=$OPTARG;;
    f)  FRAMES=$OPTARG;;
    v)  TCOUPL=$OPTARG;;
    K)  REF_T=$OPTARG;;
    b)  PCOUPL=$OPTARG;;
    Y)  PTYPE=$OPTARG;;
    B)  REF_P=$OPTARG;;
    R)  COMPRESSIBILITY=$OPTARG;;
    Z)  PBC=$OPTARG;;
    V)  SOLV_LENGTH=$OPTARG;;
    esac
done


# Calculated values based on input variables:

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"  # Directory where this script is located
INPUT_FILE=$DIR/../Structure-Files/$MONOMER  # Full path to file where monomer structure is located
Z_BOX_VECTOR=$((LAYERS+5)) # kind of arbitrary but should work
MOL_LLC=$((NO_MONOMERS*NOPORES*LAYERS))  # For topology
MOL_NA=$((NO_MONOMERS*NOPORES*LAYERS))
NSTEPS_SOLV=$(echo "$SOLV_LENGTH*1000/$DT" | bc)
NSTEPS_MD=$(echo "$SIM_LENGTH*1000/$DT" | bc)  # Number of steps to be taken during simulation to simulate the desired length of time
NSTXOUT=$(echo "$NSTEPS_MD/$FRAMES" | bc) # Information output to trajectory
NSTVOUT=$(echo "$NSTEPS_MD/$FRAMES" | bc)
NSTFOUT=$(echo "$NSTEPS_MD/$FRAMES" | bc)
NSTENERGY=$(echo "$NSTEPS_MD/$FRAMES" | bc)
NSTXOUT_SOLV=$(echo "$NSTEPS_SOLV/$FRAMES" | bc) # Information output to trajectory
NSTVOUT_SOLV=$(echo "$NSTEPS_SOLV/$FRAMES" | bc)
NSTFOUT_SOLV=$(echo "$NSTEPS_SOLV/$FRAMES" | bc)
NSTENERGY_SOLV=$(echo "$NSTEPS_SOLV/$FRAMES" | bc)
NP=$((NODES*2))

# Edit input files:

# Energy Minimization
sed -i -e "s/INTEGRATOR/${INTEGRATOR_EM}/g" em.mdp
sed -i -e "s/NSTEPS_EM/${NSTEPS_EM}/g" em.mdp
sed -i -e "s/CUTOFF_EM/${CUTOFF_EM}/g" em.mdp
sed -i -e "s/NSTLIST_EM/${NSTLIST}/g" em.mdp

# Wiggle.mdp
sed -i -e "s/SIM_TITLE/${SIM_TITLE}/g" wiggle.mdp
sed -i -e "s/CUTOFF_MD/${CUTOFF_MD}/g" wiggle.mdp
sed -i -e "s/INTEGRATOR_MD/${INTEGRATOR_MD}/g" wiggle.mdp
sed -i -e "s/DT/${DT}/g" wiggle.mdp
sed -i -e "s/NSTEPS_MD/${NSTEPS_MD}/g" wiggle.mdp
sed -i -e "s/NSTXOUT/${NSTXOUT}/g" wiggle.mdp
sed -i -e "s/NSTVOUT/${NSTVOUT}/g" wiggle.mdp
sed -i -e "s/NSTFOUT/${NSTFOUT}/g" wiggle.mdp
sed -i -e "s/NSTENERGY/${NSTENERGY}/g" wiggle.mdp
sed -i -e "s/TCOUPL/${TCOUPL}/g" wiggle.mdp
sed -i -e "s/REF_T/${REF_T}/g" wiggle.mdp
sed -i -e "s/PCOUPL/${PCOUPL}/g" wiggle.mdp
sed -i -e "s/PTYPE/${PTYPE}/g" wiggle.mdp
sed -i -e "s/REF_P/${REF_P}/g" wiggle.mdp
sed -i -e "s/COMPRESSIBILITY/${COMPRESSIBILITY}/g" wiggle.mdp
sed -i -e "s/PBC/${PBC}/g" wiggle.mdp

# Wiggle_solv.mdp
sed -i -e "s/NSTENERGY/${NSTENERGY_SOLV}/g" wiggle_solv.mdp
sed -i -e "s/SIM_TITLE/${SIM_TITLE}/g" wiggle_solv.mdp
sed -i -e "s/CUTOFF_MD/${CUTOFF_MD}/g" wiggle_solv.mdp
sed -i -e "s/INTEGRATOR_MD/${INTEGRATOR_MD}/g" wiggle_solv.mdp
sed -i -e "s/DT/${DT}/g" wiggle_solv.mdp
sed -i -e "s/NSTEPS_SOLV/${NSTEPS_SOLV}/g" wiggle_solv.mdp
sed -i -e "s/NSTENERGY/${NSTENERGY_SOLV}/g" wiggle_solv.mdp
sed -i -e "s/NSTXOUT/${NSTXOUT_SOLV}/g" wiggle_solv.mdp
sed -i -e "s/NSTVOUT/${NSTVOUT_SOLV}/g" wiggle_solv.mdp
sed -i -e "s/NSTFOUT/${NSTFOUT_SOLV}/g" wiggle_solv.mdp
sed -i -e "s/TCOUPL/${TCOUPL}/g" wiggle_solv.mdp
sed -i -e "s/REF_T/${REF_T}/g" wiggle_solv.mdp
sed -i -e "s/PCOUPL/${PCOUPL}/g" wiggle_solv.mdp
sed -i -e "s/REF_P/${REF_P}/g" wiggle_solv.mdp
sed -i -e "s/COMPRESSIBILITY/${COMPRESSIBILITY}/g" wiggle_solv.mdp

# NaPore.top
sed -i -e "s/MOL_LLC/${MOL_LLC}/g" NaPore.top
sed -i -e "s/MOL_NA/${MOL_NA}/g" NaPore.top

# NaPore_water.top
sed -i -e "s/MOL_LLC/${MOL_LLC}/g" NaPore_water.top
sed -i -e "s/MOL_NA/${MOL_NA}/g" NaPore_water.top

# Build Structure Based on user-defined inputs


python $DIR/../Structure_Builder/Structure_Builder_for_Bash.py -i $INPUT_FILE -l $LAYERS -m $NO_MONOMERS -r $RADIUS -p $PORE2PORE -n $NOPORES -d $DBWL >> initial.gro

# Put the structure in a box

gmx editconf -f initial.gro -o box.gro -c -bt triclinic -box $XVECT $YVECT $Z_BOX_VECTOR -angles 90 90 120

# Prepare input file (.tpr) for energy minimization

gmx grompp -f em.mdp -c box.gro -p NaPore.top -o box_em.tpr

# Run Energy Minimization

mpirun -np $NP gmx_mpi mdrun -v -deffnm box_em

# Extract Potential Energy from log file
ENERGY1=$(cat box_em.log | grep 'Potential Energy' | awk '{print substr($0,21,5}')

# If the potential energy is positive, then the box vector is incremented by a fixed amount until the energy comes out negative

while [ $(echo " $ENERGY1 > 0" | bc) -eq 1 ]; do
        XVECT1=$(echo "$XVECT + $INCREMENT" | bc -l)
        YVECT1=$(echo "$YVECT + $INCREMENT" | bc -l)
        gmx editconf -f initial.gro -o box.gro -c -bt triclinic -box $XVECT $YVECT $THICKNESS -angles 90 90 120
        gmx grompp -f em.mdp -c box.gro -p NaPore.top -o box_em.tpr
        mpirun -np $NP gmx_mpi mdrun -v -deffnm box_em
	ENERGY1=$(cat box_em.log | grep 'Potential Energy' | awk '{print substr($0,21,5)}')
done

# prepare input file for simulation (just letting it wiggle around)

gmx grompp -f wiggle.mdp -c box_em.gro -p NaPore.top -o wiggle_vac.tpr

# run simulation for amount of time specified in wiggle.mdp

mpirun -np $NP gmx_mpi mdrun -v -deffnm wiggle_vac

INPUT_VAC_FILE=wiggle_vac.gro

THICKNESS=$(python /$DIR/../Data-Analysis/Thickness_Bash.py -i $INPUT_VAC_FILE -w $WATER_LAYER)

# Edit box to correct dimensions (i.e. so there is room for the specified amount of water)

gmx editconf -f wiggle_vac.gro -o new_box.gro -c -bt triclinic -box 8.5 8.5 $THICKNESS -angles 90 90 60

# energy minimize in this new box since everything is shifted around

gmx grompp -f em.mdp -c new_box.gro -p NaPore.top -o em_new_box.tpr

mpirun -np $NP gmx_mpi mdrun -v -deffnm em_new_box

# Extract Potential energy from that minimization run. A positive value indicates that the box vectors are too small

ENERGY2=$(cat em_new_box.log | grep 'Potential Energy' | awk '{print substr($0,21,5)}')
# if the potential energy is positive, then the box vector is increased by a defined increment until the energy is negative after energy minimzation

while [ $(echo "  $ENERGY2 > 0" | bc) -eq 1 ]; do
        XVECT2=$(echo "$XVECT + $INCREMENT" | bc -l)
        YVECT2=$(echo "$YVECT + $INCREMENT" | bc -l)
        gmx editconf -f wiggle_vac.gro -o new_box.gro -c -bt triclinic -box $XVECT $YVECT $THICKNESS -angles 90 90 60
        gmx grompp -f em.mdp -c new_box.gro -p NaPore.top -o em_new_box.tpr
        mpirun -np $NP gmx_mpi mdrun -v -deffnm em_new_box
	ENERGY2=$(cat em_new_box.log | grep 'Potential Energy' | awk '{print substr($0,21,5)}')
done

# run simulation in new energy minimized box (Correct_box)

gmx grompp -f wiggle.mdp -c em_new_box.gro -p NaPore.top -o Correct_box.tpr

mpirun -np $NP gmx_mpi mdrun -v -deffnm Correct_box

# Now solvate the box

gmx solvate -cp Correct_box.gro -cs spc216.gro -p NaPore_water.top -o water.gro

# Now run a simulation with water in the box
# first energy minimize

gmx grompp -f em.mdp -c water.gro -p NaPore_water.top -o water_em.tpr

mpirun -np $NP gmx_mpi mdrun -v -deffnm water_em

# Now run the actual simulation

gmx grompp -f wiggle_solv.mdp -c water_em.gro -p NaPore_water.top -o water_wiggle.tpr

mpirun -np $NP gmx_mpi mdrun -v -deffnm water_wiggle

# remove unnecessary files
rm initial.gro water.gro new_box.gro #remove now unnecessary file
rm box.gro # remove unnecessary file
find . -type f -name 'Correct_box'\* -exec rm {} \;
find . -type f -name 'water_em'\* -exec rm {} \;
find . -type f -name 'em_new_box'\* -exec rm {} \;
find . -type f -name 'box_em'\* -exec rm {} \;
find . -type f -name '#'\* -exec rm {} \;
find . -type f -name 'wiggle_vac'\* -exec rm {} \;