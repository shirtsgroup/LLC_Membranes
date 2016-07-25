#!/bin/bash
# Variable Defaults:

# Approximate length of simulation: (jobs will start faster if you have a good estimate)
SIM_LENGTH_HOURS=24
SIM_LENGTH_MIN=00
SIM_LENGTH_SEC=00
NODES=4
RESOURCE='janus'
QOS='janus'  # Quality of Service
NTASKS_PER_NODE=1

#MPI Options
MPI="on"
NODES=16

# Choose which monomer to build with
MONOMER='monomer4.pdb'  # Structure file to be used
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
SIM_TITLE="'Equilibration in Vacuum'"  # Title of simulation
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
SOLVATION="off"

# Reference for flags associated with each variable

# -h  :   SIM_LENGTH_HOURS ... estimated simulation length, hours (in addition to min and seconds)
# -m  :   SIM_LENGTH_MIN ... estimated simulation length, minutes (in addition to hours and seconds)
# -s  :   SIM_LENGTH_SEC ... estimated simulation length, seconds (in addition to hours and minutes)
# -n  :   NODES ... number of nodes to use
# -q  :   QOS ... quality of service
# -N  :   NTASKS_PER_NODE ... Number of tasks assigned to each node
# -M  :   MONOMER ... Structure file to be used
# -I  :   INTEGRATOR_EM ... Integrator for energy minimization
# -S  :   NSTEPS_EM ... Maximum number of steps to take for energy minimization
# -c  :   CUTOFF_EM ... Cut-off Scheme
# -t  :   NSTLIST ... Neighborlist - changed automatically by gromacs unless it is -1- set equal to 1
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
# -Z  :   PBC ... Periodic Boundary directions-T
# -V  :   SOLV_LENGTH ... length of final simulation
# -W  :   SOLVATION ... Turn solvation on or off
# -H  :   RESOURCE ... Which machine is being using


while getopts "h:m:s:n:q:N:M:I:S:c:t:o:r:p:P:w:l:x:y:e:T:C:i:D:L:f:v:K:b:Y:B:R:Z:V:W:H:" opt; do
    case $opt in
    h)  SIM_LENGTH_HOURS=$OPTARG;;
    m)  SIM_LENGTH_MIN=$OPTARG;;
    s)  SIM_LENGTH_SEC=$OPTARG;;
    n)  NODES=$OPTARG;;
    q)  QOS=$OPTARG;;
    N)  NTASKS_PER_NODE=$OPTARG;;
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
    W)  SOLVATION=$OPTARG;;
    H)  RESOURCE=$OPTARG;;
    esac
done

if [ ${RESOURCE} == 'janus' ]; then
    NP=$((NODES*2))
    echo '#!/bin/bash' > Run_Job.sh
    echo '' >> Run_Job.sh
    echo '#SBATCH --job-name' ${SIM_TITLE} >> Run_Job.sh
    echo '#SBATCH --qos' ${QOS} >> Run_Job.sh
    echo '#SBATCH --nodes' ${NODES} >> Run_Job.sh
    echo '#SBATCH --ntasks-per-node' ${NTASKS_PER_NODE} >> Run_Job.sh
    echo '#SBATCH --time' ${SIM_LENGTH_HOURS}:${SIM_LENGTH_MIN}:${SIM_LENGTH_SEC} >> Run_Job.sh
    echo '' >> Run_Job.sh
    echo "ml slurm" >> Run_Job.sh
    echo "ml gromacs" >> Run_Job.sh
    echo "ml python/2.7.10" >> Run_Job.sh
    echo "ml numpy" >> Run_Job.sh
    echo '' >> Run_Job.sh
    echo "Build_and_Sim.sh -M ${MONOMER} -I ${INTEGRATOR_EM} -S ${NSTEPS_EM} -c ${CUTOFF_EM} -t ${NSTLIST} -o ${NO_MONOMERS} \
        -r ${RADIUS} -p ${PORE2PORE} -P ${NOPORES} -w ${DBWL} -l ${LAYERS} -x ${XVECT} -y ${YVECT} -e ${INCREMENT} \
        -T ${SIM_TITLE} -C ${CUTOFF_MD} -i ${INTEGRATOR_MD} -D ${DT} -L ${SIM_LENGTH} -f ${FRAMES} -v ${TCOUPL} -K ${REF_T}\
        -b ${PCOUPL} -Y ${PTYPE} -B ${REF_P} -R ${COMPRESSIBILITY} -Z ${PBC} -V ${SOLV_LENGTH} -n ${NODES} -s ${SOLVATION} \
        -m ${MPI}" >> Run_Job.sh
fi

if [ ${RESOURCE} == 'bridges' ]; then
    NP=$((NTASKS_PER_NODE*NODES))
    echo '#!/bin/bash' > Run_Job.sh
    echo '' >> Run_Job.sh
    echo '#SBATCH --nodes' ${NODES} >> Run_Job.sh
    echo '#SBATCH --ntasks-per-node' ${NTASKS_PER_NODE} >> Run_Job.sh
    echo '#SBATCH --time' ${SIM_LENGTH_HOURS}:${SIM_LENGTH_MIN}:${SIM_LENGTH_SEC} >> Run_Job.sh
    echo '' >> Run_Job.sh
    echo "module load gromacs" >> Run_Job.sh
    echo "module load python/2.7.11_gcc" >> Run_Job.sh
    echo "source /home/bjc/Programs/Gromacs/bin/GMXRC" >> Run_Job.sh
    echo '' >> Run_Job.sh
    echo "Build_and_Sim.sh -M ${MONOMER} -I ${INTEGRATOR_EM} -S ${NSTEPS_EM} -c ${CUTOFF_EM} -t ${NSTLIST} -o ${NO_MONOMERS} \
        -r ${RADIUS} -p ${PORE2PORE} -P ${NOPORES} -w ${DBWL} -l ${LAYERS} -x ${XVECT} -y ${YVECT} -e ${INCREMENT} \
        -T ${SIM_TITLE} -C ${CUTOFF_MD} -i ${INTEGRATOR_MD} -D ${DT} -L ${SIM_LENGTH} -f ${FRAMES} -v ${TCOUPL} -K ${REF_T}\
        -b ${PCOUPL} -Y ${PTYPE} -B ${REF_P} -R ${COMPRESSIBILITY} -Z ${PBC} -V ${SOLV_LENGTH} -n ${NODES} -s ${SOLVATION} \
        -m ${MPI} -H ${RESOURCE} -a ${NTASKS_PER_NODE}" >> Run_Job.sh
else
    echo "Specify a valid resource"
    echo "NOTE: names are all lowercase"
    echo "i.e. don't type Janus or Bridges or you might get this message"
fi
