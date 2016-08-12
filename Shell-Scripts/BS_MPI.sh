#!/bin/bash
# Script to build a structure to a specified number of layers using a specified
# monomer, then energy minimize the structure, and finally run a simulation
# based on inputs to the .mdp files

# Copy in necessary files
ifv.sh

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
SIM_TITLE='Equilibration in Vacuum'  # Title of simulation
CUTOFF_MD='verlet'  # Cut-off scheme for simulation
INTEGRATOR_MD='md'  # Integrator type for simulation
DT=0.002  # Time step (ps)
SIM_LENGTH=1  # nanoseconds
FRAMES=50  # Number of frames in trajectory
TCOUPL='v-rescale'
REF_T=300  # Reference Temperature, K
PCOUPL='berendsen'
PTYPE='semiisotropic'
REF_P=1  # Reference Pressure, bar
COMPRESSIBILITY=4.5e-5  # Isothermal compressibility, bar^-1
PBC='xyz'

# Reference for flags associated with each variable

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
# -L  :   SIM_LENGTH ... Simulation length, nanoseconds
# -f  :   FRAMES ... Number of frames in trajectory
# -v  :   TCOUPL ... Temperature Coupling
# -K  :   REF_T ... Reference Temperature, K
# -b  :   PCOUPL ... Pressure Coupling
# -Y  :   PTYPE ... i.e. Isotropic, semiisotropic
# -B  :   REF_P ... Reference Pressure, bar
# -R  :   COMPRESSIBILITY ... Isothermal compressibility, bar^-1
# -Z  :   PBC ... Periodic Boundary directions
# -m  :   NODES ... Number of nodes to use on parallel machine

while getopts "M:I:S:c:t:o:r:p:P:w:l:x:y:e:T:C:i:D:L:f:v:K:b:Y:B:R:Z:m:" opt; do
    case $opt in
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
    m)  NODES=$OPTARG;;
    esac
done

# Calculated values based on input variables:

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"  # Directory where this script is located
#INPUT_FILE=$DIR/../Structure-Files/$MONOMER  # Full path to file where monomer structure is located
Z_BOX_VECTOR=$((LAYERS+5)) # kind of arbitrary but should work
MOL_LLC=$((NO_MONOMERS*NOPORES*LAYERS))  # For topology
MOL_NA=$((NO_MONOMERS*NOPORES*LAYERS))
NSTEPS_MD=$(echo "$SIM_LENGTH*1000/$DT" | bc)  # Number of steps to be taken during simulation to simulate the desired length of time
NSTXOUT=$(echo "$NSTEPS_MD/$FRAMES" | bc) # Information output to trajectory
NSTVOUT=$(echo "$NSTEPS_MD/$FRAMES" | bc)
NSTFOUT=$(echo "$NSTEPS_MD/$FRAMES" | bc)
NSTENERGY=$(echo "$NSTEPS_MD/$FRAMES" | bc)
NP=$((NODES*2))  # for janus ... needs to be modified for other machines once parallelization scheme is known

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

# Edit Topology
sed -i -e "s/MOL_LLC/${MOL_LLC}/g" NaPore.top
sed -i -e "s/MOL_NA/${MOL_NA}/g" NaPore.top

# Build Structure Based on user-defined inputs

python $DIR/../Structure_Builder/Orient_Plane.py -i $MONOMER -l $LAYERS -m $NO_MONOMERS -r $RADIUS -p $PORE2PORE -n $NOPORES -d $DBWL >> initial.gro

# Put the structure in a box

gmx editconf -f initial.gro -o box.gro -c -bt triclinic -box $XVECT $YVECT $Z_BOX_VECTOR -angles 90 90 120

# Prepare input file (.tpr) for energy minimization

gmx grompp -f em.mdp -c box.gro -p NaPore.top -o box_em.tpr

# Run Energy Minimization

mpirun -np $NP gmx_mpi mdrun -v -deffnm box_em

# Extract Potential Energy from log file
ENERGY=$(cat box_em.log | grep 'Potential Energy' | awk '{print substr($0,21,5}')

# If the potential energy is positive, then the box vector is incremented by a fixed amount until the energy comes out negative

while [ $(echo " $ENERGY > 0" | bc) -eq 1 ]; do
        XVECT=$(echo "$XVECT + $INCREMENT" | bc -l)
        YVECT=$(echo "$YVECT + $INCREMENT" | bc -l)
        gmx editconf -f initial.gro -o box.gro -c -bt triclinic -box $XVECT $YVECT $Z_BOX_VECTOR -angles 90 90 120
        gmx grompp -f em.mdp -c box.gro -p NaPore.top -o box_em.tpr
        mpirun -np $NP gmx_mpi mdrun -v -deffnm box_em
        ENERGY=$(cat box_em.log | grep 'Potential Energy' | awk '{print substr($0,21,5)}')
done

# prepare input file for simulation (just letting it wiggle around)

gmx grompp -f wiggle.mdp -c box_em.gro -p NaPore.top -o wiggle.tpr

# run simulation for amount of time specified in wiggle.mdp

mpirun -np $NP gmx_mpi mdrun -v -deffnm wiggle

gmx trjconv -f wiggle.trr -s wiggle.tpr -o traj.gro
# remove unecessary files
rm box_em*
rm \#*
rm initial.gro
rm box.gro


