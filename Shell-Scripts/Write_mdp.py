#!/usr/bin/python

import os
import argparse

parser = argparse.ArgumentParser(description = 'Write .mdp files for a vacuum simulation')
# Flags mirror those input into bash script for running simulation
parser.add_argument('-T', '--title', default='Vacuum Simulation', help = 'Title for the simulation')
parser.add_argument('-C', '--CUTOFF_MD', default='verlet', help = 'Cut-off algorithm')
parser.add_argument('-i', '--INTEGRATOR_MD', default='md')
parser.add_argument('-D', '--DT', default=0.002, help = 'Time step (picoseconds)')
parser.add_argument('-v', '--TCOUPL', default='v-rescale', help = 'Temperature Coupling Algorithm')
parser.add_argument('-K', '--REF_T', default=300, help = 'Reference Temperature')
parser.add_argument('-b', '--PCOUPL', default='berendsen', help = 'Pressure Coupling Algorithm')
parser.add_argument('-Y', '--PTYPE', default='semiisotropic', help = 'Type of pressure coupling')
parser.add_argument('-B', '--REF_P', default=1, help = 'Reference Pressure (bar)')
parser.add_argument('-R', '--COMPRESSIBILITY', default='4.5e-5', help = 'Isothermal Compressibility, bar^-1')
parser.add_argument('-Z', '--PBC', default='xyz', help = 'Periodic Boundary Condition direction')
parser.add_argument('-L', '--SIM_LENGTH', default=1, help = 'Length of Vacuum Simulation (ns)')
parser.add_argument('-I', '--INTEGRATOR_EM', default = 'steep', help= 'Integrator used for energy minimization')
parser.add_argument('-S', '--NSTEPS_EM', default = 50000, help = 'Maximum number of steps to take during energy minimization')
parser.add_argument('-c', '--CUTOFF_EM', default = 'verlet', help = 'Cutoff algorithm for energy minimzation')
parser.add_argument('-t', '--NSTLIST', default = 40, help = 'Neighbor search list')
parser.add_argument('-f', '--FRAMES', default = 50, help = 'Number of Frames')
parser.add_argument('-m', '--monomer', default = 'LLC', help = 'Monomer which structure is built with')
parser.add_argument('-s', '--solvated', default = 'off', help = 'Will create a .mdp file for solvated system if this is set to "on"')
parser.add_argument('-V', '--SOLV_LENGTH', default = 1, help = 'Length of simulation of solvated system')


args = parser.parse_args()

monomer = args.monomer

if monomer == 'LLC':
    grps = ['LLC', 'NA']
    grps_solv = ['LLC', 'NA', 'SOL']
elif monomer == 'BCC':
    grps = ['BCC', 'BR']
    grps_solv = ['BCC', 'BR', 'SOL']

# Energy minimization .mdp file
title = 'title = Energy Minimization'
integrator = 'integrator = %s' %args.INTEGRATOR_EM
nsteps = 'nsteps = %s' %args.NSTEPS_EM
cutoff_scheme = 'cutoff-scheme = %s' %args.CUTOFF_EM
nstlist = 'nstlist = %s' %args.NSTLIST

f1 = open('em.mdp', 'w')
f1.writelines([title + '\n', integrator + '\n', nsteps + '\n', cutoff_scheme + '\n', nstlist + '\n'])

# Vacuum Simulation .mdp file
title = 'title = %s' %args.title
cutoff_scheme = 'cutoff-scheme = %s' %args.CUTOFF_MD
integrator = 'integrator = %s' %args.INTEGRATOR_MD
dt = 'dt = %s' %args.DT
steps = int(float(args.SIM_LENGTH) * 1000 / float(args.DT))
nsteps = 'nsteps = %s' %steps
nstxout = 'nstxout = %s' %(int(steps/int(args.FRAMES)))
nstvout = 'nstvout = %s' %(int(steps/int(args.FRAMES)))
nstfout = 'nstfout = %s' %(int(steps/int(args.FRAMES)))
nstenergy = 'nstenergy = %s' %(int(steps/int(args.FRAMES)))
nstlist = 'nstlist = %s' %args.NSTLIST
tcoupl = 'Tcoupl = %s' %args.TCOUPL
tc_grps = 'tc_grps = %s' %' '.join(grps)
tau_t = 'tau_t = %s' %' '.join([str(0.1) for i in grps])
ref_t = 'ref_t = %s' %' '.join([str(args.REF_T) for i in grps])
Pcoupl = 'Pcoupl = %s' %args.PCOUPL
Pcoupltype = 'Pcoupltype = %s' %args.PTYPE
if args.PTYPE == 'semiisotropic':
    compress = 'compressibility = %s' %' '.join([str(args.COMPRESSIBILITY), '0'])
    ref_p = 'ref_p = %s' %' '.join([str(args.REF_P) for i in grps])
else:
    compress = 'compressibility = %s' %args.COMPRESSIBILITY
    ref_p = 'ref_p = %s' %args.REF_P
pbc = 'pbc = %s' %args.PBC

f2 = open('wiggle.mdp', 'w')
f2.writelines([title + '\n', cutoff_scheme + '\n', integrator + '\n', dt + '\n', nsteps + '\n', 'continuation = no\n',
               'constraints = all-bonds\n', 'constraint-algorithm = lincs\n', 'lincs-iter = 1\n', 'lincs-order = 4\n',
               nstxout + '\n', nstvout + '\n', nstfout + '\n', nstenergy + '\n', nstlist + '\n', 'ns_type = grid\n'
               'rlist = 1.2\n', 'rcoulomb = 1.2\n', 'rvdw = 1.2\n', 'coulombtype = PME\n', 'pme_order = 4\n',
               'fourierspacing = 0.16\n', tcoupl + '\n', tc_grps + '\n', tau_t + '\n', ref_t + '\n', Pcoupl + '\n', Pcoupltype + '\n',
               ref_p + '\n', compress + '\n', 'gen_vel = no\n', pbc + '\n', 'DispCorr = Ener\n'])

if args.solvated == 'on':
    title = 'title = Solvated System'
    steps = int(float(args.SOLV_LENGTH) * 1000 / float(args.DT))
    nsteps = 'nsteps = %s' %steps
    nstxout = 'nstxout = %s' %(steps/int(args.FRAMES))
    nstvout = 'nstvout = %s' %(steps/int(args.FRAMES))
    nstfout = 'nstfout = %s' %(steps/int(args.FRAMES))
    nstenergy = 'nstenergy = %s' %(steps/int(args.FRAMES))
    tcoupl = 'Tcoupl = %s' %args.TCOUPL
    tc_grps = 'tc_grps = %s' %' '.join(grps_solv)
    tau_t = 'tau_t = %s' %' '.join([str(0.1) for i in grps_solv])
    ref_t = 'ref_t = %s' %' '.join([str(args.REF_T) for i in grps_solv])
    Pcoupl = 'Pcoupl = %s' %args.PCOUPL
    Pcoupltype = 'Pcoupltype = isotropic'
    compress = 'compressibility = %s' %args.COMPRESSIBILITY
    ref_p = 'ref_p = %s' %args.REF_P
    pbc = 'pbc = %s' %args.PBC

    f3 = open('wiggle_solv.mdp', 'w')
    f3.writelines([title + '\n', cutoff_scheme + '\n', integrator + '\n', dt + '\n', nsteps + '\n', 'continuation = no\n',
                   'constraints = all-bonds\n', 'constraint-algorithm = lincs\n', 'lincs-iter = 1\n', 'lincs-order = 4\n',
                   nstxout + '\n', nstvout + '\n', nstfout + '\n', nstenergy + '\n', nstlist + '\n', 'ns_type = grid\n'
                   'rlist = 1.2\n', 'rcoulomb = 1.2\n', 'rvdw = 1.2\n', 'coulombtype = PME\n', 'pme_order = 4\n',
                   'fourierspacing = 0.16\n', tcoupl + '\n', tc_grps + '\n', tau_t + '\n', ref_t + '\n', Pcoupl + '\n', Pcoupltype + '\n',
                   ref_p + '\n', compress + '\n', 'gen_vel = no\n', pbc + '\n', 'DispCorr = Ener\n'])
