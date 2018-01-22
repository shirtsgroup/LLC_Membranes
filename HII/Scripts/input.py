#! /usr/bin/env python

"""
Write input files. Default settings are not included. Add an argument if you need to change a default that isn't shown
"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
import os
import argparse
import lc_class
import bcc_class


def initialize():

    parser = argparse.ArgumentParser(description='Write .mdp files and topology files for various types of simulation')

    parser.add_argument('-T', '--title', default='Generic Molecular Dynamics Simulation', type=str, help='Simulation Title')
    parser.add_argument('-b', '--build_mon', default='NAcarb11V', type=str, help='Monomer structure used for build')
    parser.add_argument('-t', '--itp', default='dipole.itp', type=str, help='Name of .itp describing monomers')
    parser.add_argument('-s', '--em_steps', default=5000, type=int, help='Steps to take during energy minimization')
    parser.add_argument('-e', '--ensemble', default='npt', type=str, help='Thermodynamic ensemble to put system in')
    parser.add_argument('-d', '--dt', default=0.002, type=float, help='time step (ps)')
    parser.add_argument('-l', '--length', default=1000, type=float, help='simulation length (ps)')
    parser.add_argument('-f', '--frames', default=500, type=int, help='number of frames')
    parser.add_argument('-p', '--pcoupltype', default='semiisotropic', type=str, help='Pressure Couple Type')
    parser.add_argument('--restraints', help='If restraints are on, another mdp option needs to be turned on, so specify '
                                             'this flag', action="store_true")
    parser.add_argument('-x', '--xlink', help='Turn this to "on" if the the system is crosslinked', action="store_true")
    parser.add_argument('-c', '--coord', default='initial.gro', type=str, help='coordinate file of system to be simulated')
    parser.add_argument('-S', '--solvate', help='Specify this if the system has water so an extra line can be added to the '
                                                'topology', action="store_true")
    parser.add_argument('--temp', default=300, help='Specify temperature at which to run simulation')
    parser.add_argument('--mdp', action="store_true", help='Only the .mdp will be written if this option is specified')
    parser.add_argument('--barostat', default='berendsen', type=str, help='pressure coupling scheme to use')
    parser.add_argument('--genvel', default='yes', type=str, help='generate velocities according to a maxwell'
                                                                     'distribution')
    parser.add_argument('--bcc', action="store_true", help='Generate input files using bicontinuous cubic files')  # probably best to reorganize the repository
    parser.add_argument('--solvent', default='water', help='Name of solvent')

    args = parser.parse_args()

    return args

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in

if __name__ == "__main__":

    args = initialize()

    if args.bcc:
        props = bcc_class.LC('%s.gro' % args.build_mon)
    else:
        props = lc_class.LC('%s.gro' % args.build_mon)

    mon_name = props.residues[0]

    if args.xlink:
        mon_top = '#include "%s/crosslinked_new.itp"' % os.getcwd()
    else:
        if args.bcc:
            mon_top = '#include "%s/../../BCC/top/topologies/%s.itp' % (location, args.build_mon)
        else:
            mon_top = '#include "%s/../top/Monomer_Tops/%s.itp' % (location, args.build_mon)

    gaff = '#include "%s/../top/Forcefields/gaff' % location  # generalized amber force field

    # Energy minimization .mdp file
    title = 'title = Energy Minimization'
    integrator = 'integrator = steep'
    nsteps = 'nsteps = %s' % args.em_steps
    cutoff_scheme = 'cutoff-scheme = verlet'
    nstlist = 'nstlist = 40'

    f = open('em.mdp', 'w')
    f.writelines([title + '\n', integrator + '\n', nsteps + '\n', cutoff_scheme + '\n', nstlist + '\n'])

    if args.ensemble == 'npt':

        a = []
        a.append(['title = NPT simulation of %s at %s K\n' % (args.build_mon, args.temp)])
        # a.append(['cutoff-scheme = verlet'])  # I think verlet is default
        a.append(['integrator = md\n'])  # this also might be default
        a.append(['dt = %s\n' % args.dt])
        a.append(['nsteps = %s\n' % int(args.length / args.dt)])
        a.append(['continuation = no\n'])
        a.append(['constraints = h-bonds\n'])
        a.append(['constraint-algorithm = lincs\n'])
        a.append(['cutoff-scheme = Verlet\n'])
        a.append(['nstxout = %s\n' % int(args.length / (args.dt * args.frames))])
        a.append(['nstvout = %s\n' % int(args.length / (args.dt * args.frames))])
        a.append(['nstfout = %s\n' % int(args.length / (args.dt * args.frames))])
        a.append(['nstenergy = %s\n' % int(args.length / (args.dt * args.frames))])
        a.append(['nstlist = 40\n'])
        a.append(['nstype = grid\n'])
        a.append(['vdwtype = PME\n'])
        a.append(['coulombtype = PME\n'])
        a.append(['Tcoupl = v-rescale\n'])
        a.append(['tc-grps = system\n'])
        a.append(['tau-t = %s\n' % str(0.1)])
        a.append(['ref-t = %s\n' % args.temp])
        if not args.restraints:
            a.append(['Pcoupl = %s\n' % args.barostat])
            a.append(['Pcoupltype = %s\n' % args.pcoupltype])
            if args.barostat == 'Parrinello-Rahman':
                a.append(['tau-p = 10\n'])  # tau-p should be at least  20 times larger than nstpcouple*dt. nstpcoupl defaults to the value of nstlist (40)
            if args.pcoupltype == 'Isotropic':
                a.append(['ref-p = 1\n'])
                a.append(['compressibility = 4.5e-5\n'])
            else:
                a.append(['ref-p = %s\n' % ' '.join([str(1) for i in props.residues])])
                a.append(['compressibility = 4.5e-5 4.5e-5\n'])
        if args.genvel == 'yes':
            a.append(['gen-vel = yes\n'])
            a.append(['gen-temp = %s\n' % args.temp])
        else:
            a.append(['gen-vel = no\n'])
        a.append(['pbc = xyz\n'])
        a.append(['DispCorr = EnerPres\n'])
        if args.xlink:
            a.append('periodic-molecules = yes\n')
            a.append('lincs-iter=4')
        if args.restraints:
            a.append(['refcoord-scaling = all\n'])

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
        a.append(['vdwtype = PME\n'])
        a.append(['coulombtype = PME\n'])
        a.append(['Tcoupl = v-rescale\n'])
        a.append(['tc_grps = system\n'])
        a.append(['tau_t = %s\n' % str(0.1)])
        a.append(['ref_t = %s\n' % args.temp])
        a.append(['gen_vel = yes\n'])
        a.append(['pbc = xyz\n'])
        a.append(['DispCorr = EnerPres\n'])
        if args.xlink:
            a.append('periodic-molecules = yes\n')
        if args.restraints:
            a.append(['refcoord_scaling = all\n'])

        f = open('%s.mdp' % args.ensemble, 'w')
        for line in a:
            f.write(line[0])

        f.close()

    if args.mdp:
        exit()

    # Create topology
    a = []
    a.append(';Forcefield\n')
    a.append('%s/gaff.itp"\n' % gaff)
    a.append('\n')
    a.append(';Monomer Topology\n')
    if args.restraints:
        a.append('#include "./%s"\n' % args.itp)
    else:
        a.append('%s"\n' % mon_top)
    a.append('\n')
    if props.no_ions > 0:
        a.append(';Ion Topology\n')
        a.append('%s/ions.itp"\n' % gaff)
        a.append('\n')
    if args.solvate:
        if args.solvent == 'water':
            a.append(';Water Topology\n')
            a.append('%s/tip3p.itp\n' % gaff)
            a.append('\n')
        elif args.solvent == 'glycerol':
            a.append(';Glycerol Topology\n')
            a.append('#include "%s/../../BCC/top/topologies/glycerol.itp\n' %location)
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
    for i in range(2, len(gro)):
        if gro[i].count('%s' % mon_name) != 0:
            nres += 1

    if props.ions:
        nion = 0
        for i in range(2, len(gro)):
            if gro[i].count('%s' % props.ions[0]) != 0:
                nion += 1

    if args.bcc:
        nmon = int(nres / (props.natoms - props.no_ions))
    else:
        nmon = int(nres / props.natoms)

    if args.restraints:
        a.append('%s                1\n' % mon_name)
    else:
        if args.xlink:
            nmon = 1
        a.append('%s                 %s\n' % (mon_name, nmon))

    if props.ions:
        a.append('%s                 %s\n' % (props.ions[0], nion))

    if args.solvate:
        if not args.bcc:
            sol = 0
            for i in range(2, len(gro)):
                if gro[i].count('SOL') != 0:
                    sol += 1

            a.append('SOL                %s\n' % int(sol/3))

    f = open('topol.top', 'w')
    for line in a:
        f.write(line)

    f.close()