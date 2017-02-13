#!/usr/bin/python

"""
Write input files. Default settings are not included. Add an argument if you need to change a default that isn't shown
"""

import os
import argparse
import LC_class

parser = argparse.ArgumentParser(description='Write .mdp files and topology files for various types of simulation')

parser.add_argument('-T', '--title', default='Generic Molecular Dynamics Simulation', type=str, help='Simulation Title')
parser.add_argument('-b', '--build_mon', default='NAcarb11Vd', type=str, help='Monomer structure used for build')
parser.add_argument('-s', '--em_steps', default=50000, type=int, help='Steps to take during energy minimization')
parser.add_argument('-e', '--ensemble', default='npt', type=str, help='Thermodynamic ensemble to put system in')
parser.add_argument('-d', '--dt', default=0.002, type=float, help='time step (ps)')
parser.add_argument('-l', '--length', default=1000, type=int, help='simulation length (ps)')
parser.add_argument('-f', '--frames', default=50, type=int, help='number of frames')
parser.add_argument('-p', '--pcoupltype', default='semiisotropic', type=str, help='Pressure Couple Type')
parser.add_argument('-r', '--restraints', default='off', type=str, help='If restraints are on, another mdp option needs'
                                                                        'to be turned on')
parser.add_argument('-x', '--xlink', default='off', type=str, help='Turn this to "on" if the the system is crosslinked')
parser.add_argument('-c', '--coord', default='initial.gro', type=str, help='coordinate file of system to be simulated')

args = parser.parse_args()

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in

exec "build_mon = LC_class.%s.build_mon" % args.build_mon
exec "ion = LC_class.%s.counterion" % args.build_mon
exec "mon_name = LC_class.%s.name" % args.build_mon
exec "no_ions = LC_class.%s.valence" % args.build_mon
exec "name = LC_class.%s.name" % args.build_mon
exec "grps = LC_class.%s.residues" % args.build_mon
exec "natoms = LC_class.%s.atoms" % args.build_mon

if args.xlink == 'on':
    mon_top = '#include "%s/crosslinked_new.itp"' % os.getcwd()
else:
    mon_top = '#include "%s/../Structure-Files/Monomer_Tops/%s.itp' % (location, args.build_mon)

gaff = '#include "%s/../Structure-Files/Forcefields/gaff/' % location  # generalized amber force field

# Energy minimization .mdp file
title = 'title = Energy Minimization'
integrator = 'integrator = steep'
nsteps = 'nsteps = %s' %args.em_steps
cutoff_scheme = 'cutoff-scheme = verlet'
nstlist = 'nstlist = 40'

f = open('em.mdp', 'w')
f.writelines([title + '\n', integrator + '\n', nsteps + '\n', cutoff_scheme + '\n', nstlist + '\n'])

if args.ensemble == 'npt':

    a = []
    a.append(['title = NPT simulation of %s\n' % args.build_mon])
    # a.append(['cutoff-scheme = verlet'])  # I think verlet is default
    a.append(['integrator = md\n'])  # this also might be default
    a.append(['dt = %s\n' % args.dt])
    a.append(['nsteps = %s\n' % int(args.length / args.dt)])
    a.append(['continuation = no\n'])
    a.append(['constraints = h-bonds\n'])
    a.append(['constraint-algorithm = lincs\n'])
    a.append(['nstxout = %s\n' % int(args.length / (args.dt * args.frames))])
    a.append(['nstvout = %s\n' % int(args.length / (args.dt * args.frames))])
    a.append(['nstfout = %s\n' % int(args.length / (args.dt * args.frames))])
    a.append(['nstenergy = %s\n' % int(args.length / (args.dt * args.frames))])
    a.append(['nstlist = 40\n'])
    a.append(['nstype = grid\n'])
    a.append(['coulombtype = PME\n'])
    a.append(['Tcoupl = v-rescale\n'])
    a.append(['tc_grps = %s\n' % ' '.join(grps)])
    a.append(['tau_t = %s\n' % ' '.join([str(0.1) for i in grps])])
    a.append(['ref_t = %s\n' % ' '.join([str(300) for i in grps])])
    a.append(['Pcoupl = berendsen\n'])
    a.append(['Pcoupltype = %s\n' % args.pcoupltype])
    a.append(['ref_p = %s\n' % ' '.join([str(1) for i in grps])])
    if args.pcoupltype == 'Isotropic':
        a.append(['compressibility = 4.5e-5\n'])
    else:
        a.append(['compressibility = 4.5e-5 0\n'])
    a.append(['gen_vel = yes\n'])
    a.append(['pbc = xyz\n'])
    a.append(['DispCorr = Ener\n'])
    if args.xlink == 'on':
        a.append('periodic-molecules = yes\n')
    if args.restraints == 'on':
        a.append(['refcoord_scaling = all\n'])

    f = open('%s.mdp' % args.ensemble, 'w')
    for line in a:
        f.write(line[0])

    f.close()

if args.ensemble == 'nvt':

    a = []
    a.append(['title = NVT simulation of %s\n' % args.build_mon])
    # a.append(['cutoff-scheme = verlet'])  # I think verlet is default
    a.append(['integrator = md\n'])  # this also might be default
    a.append(['dt = %s\n' % args.dt])
    a.append(['nsteps = %s\n' % int(args.length / args.dt)])
    a.append(['continuation = no\n'])
    a.append(['constraints = h-bonds\n'])
    a.append(['constraint-algorithm = lincs\n'])
    a.append(['nstxout = %s\n' % int(args.length / (args.dt * args.frames))])
    a.append(['nstvout = %s\n' % int(args.length / (args.dt * args.frames))])
    a.append(['nstfout = %s\n' % int(args.length / (args.dt * args.frames))])
    a.append(['nstenergy = %s\n' % int(args.length / (args.dt * args.frames))])
    a.append(['nstlist = 40\n'])
    a.append(['nstype = grid\n'])
    a.append(['coulombtype = PME\n'])
    a.append(['Tcoupl = v-rescale\n'])
    a.append(['tc_grps = %s\n' % ' '.join(grps)])
    a.append(['tau_t = %s\n' % ' '.join([str(0.1) for i in grps])])
    a.append(['ref_t = %s\n' % ' '.join([str(300) for i in grps])])
    a.append(['gen_vel = yes\n'])
    a.append(['pbc = xyz\n'])
    a.append(['DispCorr = Ener\n'])
    if args.xlink == 'on':
        a.append('periodic-molecules = yes\n')
    if args.restraints == 'on':
        a.append(['refcoord_scaling = all\n'])

    f = open('%s.mdp' % args.ensemble, 'w')
    for line in a:
        f.write(line[0])

    f.close()

a = []
a.append(';Forcefield\n')
a.append('%s\n' % gaff)
a.append('\n')
a.append(';Monomer Topology\n')
if args.restraints == 'on':
    a.append('#include "./dipole.itp"\n')
else:
    a.append('%s/gaff.itp\n' % mon_top)
a.append('\n')
a.append(';Ion Topology\n')
a.append('%s/ions.itp\n' % gaff)
a.append('\n')
a.append('[ system ]\n')
a.append('%s simulation of %s\n' % (args.ensemble, args.build_mon))
a.append('\n')
a.append('[ molecules ]\n')
a.append('; Compound         nmols\n')

f = open('%s' % args.coord)
gro = []
for line in f:
    gro.append(line)

nres = 0
nion = 0
for i in range(len(gro)):
    if gro[i].count('%s' % mon_name) != 0:
        nres += 1
    if gro[i].count('%s' % ion) != 0:
        nion += 1

nmon = int(nres / natoms)

if args.restraints == 'off':
    a.append('%s                %s\n' % (mon_name, nmon))
else:
    a.append('%s                1\n' % (mon_name))
a.append('%s                 %s\n' % (ion, nion))

f = open('topol.top', 'w')
for line in a:
    f.write(line)

f.close()