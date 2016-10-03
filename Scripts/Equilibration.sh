#!/bin/bash

# Procedure for equilibration, especially for tricky systems

TOP="NaPore.top"
LENGTH1=0.01 # ns
TAU_T=0.5
TAU_P=10
STEP=0.0005
EMSTEPS=5000

while getopts "g:t:l:u:e:" opt; do
    case ${opt} in
    g) GRO=$OPTARG;;
    t) TOP=$OPTARG;;
    l) LENGTH1=$OPTARG;;
    u) TAU_T=$OPTARG;;
    s) STEP=$OPTARG;;
    e) EMSTEPS=$OPTARG;;
    esac
done

# Write input files : wiggle.mdp, em.mdp, em_lbfgs.mdp, NaPore.top
Write_Input.py -D ${STEP} -L ${LENGTH1} -e "NVT" -u ${TAU_T} -S ${EMSTEPS} -g 'yes' -a 6 -A 8 -d 'alternating'

gmx grompp -f em.mdp -p ${TOP} -c ${GRO} -o em

gmx mdrun -v -deffnm em  # run steepest descent energy minimization

gmx grompp -f em_lbfgs.mdp -p ${TOP} -c ${GRO} -o em_lbfgs

gmx mdrun -v -deffnm em_lbfgs  # continue energy minimization where it left off using l-bfgs integrator

gmx grompp -f wiggle.mdp -p ${TOP} -c em_lbfgs.gro -o wiggle

gmx mdrun -v -deffnm wiggle  # Run a short simulation

echo "13" | gmx energy -f wiggle.edr  # Get Temperature vs. time data

nohup Plot_xvg.py &  # Plot T vs. Time. Detach the graph completely using the nohup/& combo

PID=$!  # The process ID of the graph that is shown

echo "Is tau_t sufficient based on this graph"
select yn in "Yes" "No"; do
    case $yn in
        Yes) kill $PID; break;;
        No ) kill $PID
             TAU_T=$(echo "scale=3; ${TAU_T}/2" | bc -l) # cuts tau_t in half
             line=$(awk '/tau_t/{ print NR; exit }' wiggle.mdp)  # gets line number where tau_t is specified
             new_line="tau_t = ${TAU_T} ${TAU_T}"  # prepare a new line to replace the old one with new tau_t's
             sed -i "${line}s/.*/${new_line}/" wiggle.mdp  # now use the new line to replace the old line
             gmx grompp -f wiggle.mdp -p ${TOP} -c em_lbfgs.gro -o wiggle  # now rerun the simulation
             gmx mdrun -v -deffnm wiggle
             echo "How about now?"
             nohup Plot_xvg.py &
             PID=$!;;
    esac
done

# Run NVT at dt = 0.002
dt_line=$(awk '/dt/{ print NR; exit }' wiggle.mdp)  # find line where dt is specified in wiggle.mdp
sed -i "${dt_line}s/.*/dt = 0.002/" wiggle.mdp  # replace line with new dt

# Run the NVT simulation for 50 ps. We need to adjust the number of steps now
new_steps=$(echo "50/${dt}" | bc)
steps_line=$(awk '/nsteps/{ print NR; exit }' wiggle.mdp)
sed -i "${steps_line}s/.*/nsteps = ${new_steps}/" wiggle.mdp

# Make atom level input file
gmx grompp -f wiggle.mdp -p NaPore.top -c wiggle.gro -o wiggle_NVT

# Run NVT simulation
gmx mdrun -v -deffnm wiggle_NVT

echo "Check for structural issues. If modifications are made, replace the .gro file with a .gro file of the same name.
       You may need to perform an energy minimization first"
echo "Continue equilibration?"
echo "1) yes = continue"
echo "2) no  = abort equilibration"
select yn in "Yes" "No"; do
    case $yn in
        Yes) break;;
        No ) exit;;
    esac
done

# Now prepare for an NPT simulation with a small time step, short tau and relatively long tau_p (10 - 50 ps)
# Add lines to NVT .mdp file so that it becomes NPT
echo "Pcoupl = berendsen" >> wiggle.mdp
echo "Pcoupltype = semiisotropic" >> wiggle.mdp
echo "ref_p = 1 1" >> wiggle.mdp
echo "compressibility = 4.5e-5 0" >> wiggle.mdp

# Assemble atomic level input file
gmx grompp -f wiggle.mdp -p NaPore.top -c wiggle_NVT.gro -o wiggle_NPT

# Run NPT simulation
gmx mdrun -v -deffnm wiggle_NPT

# Now run a normal NPT simulation (i.e. change the time step back to 0.002