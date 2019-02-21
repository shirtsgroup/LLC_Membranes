#!/usr/bin/env python

import argparse
import os
import subprocess
from LLC_Membranes.setup.gentop import SystemTopology
from LLC_Membranes.setup.genmdp import SimulationMdp

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in


def initialize():

    parser = argparse.ArgumentParser(description='Solvate system and adjust so there is a certain wt % of water')

    parser.add_argument('-i', '--initial', default='initial.gro', help='Coordinate file to add water to')
    parser.add_argument('-ratio', '--ratio', default=1.5, type=float, help='Ratio of water in pores to water in tails')
    parser.add_argument('-wt', '--weight_percent', default=10, type=float, help='Total weight percent of water')
    parser.add_argument('-tol', '--tolerance', default=1, type=int, help='Number of water molecules')
    parser.add_argument('-guess_range', default=[.4, 1], nargs='+', help='If water_content.db has no entries for the '
                        'build monomer, an initial radius will be randomly selected from this range')
    parser.add_argument('-guess_stride', default=0.2, type=float, help='How far above/below the highest/lowest value to'
                        'make the next guess at pore radius if you need more/less water than the bounds of the water '
                        'content database (nm)')
    parser.add_argument('-o', '--output', default='solvated_final.gro', help='Name of output file')

    # parallelization
    parser.add_argument('-mpi', '--mpi', action="store_true", help='Run MD simulations in parallel')
    parser.add_argument('-np', '--nproc', default=4, help='Number of MPI processes')

    # same flags as to build.py
    parser.add_argument('-b', '--build_monomer', default='NAcarb11V.gro', type=str, help='Name of single monomer'
                        'structure file (.gro format) used to build full system')
    parser.add_argument('-m', '--monomers_per_column', default=20, type=int, help='Number of monomers to stack in each'
                                                                                  'column')
    parser.add_argument('-c', '--ncolumns', default=5, type=int, help='Number of columns used to build each pore')
    parser.add_argument('-r', '--pore_radius', default=.6, type=float, help='Initial guess at pore radius (nm)')
    parser.add_argument('-p', '--p2p', default=4.5, type=float, help='Initial pore-to-pore distance (nm)')
    parser.add_argument('-n', '--nopores', default=4, type=int, help='Number of pores (only works with 4 currently)')
    parser.add_argument('-d', '--dbwl', default=.37, type=float, help='Distance between vertically stacked monomers'
                                                                      '(nm)')
    parser.add_argument('-pd', '--parallel_displaced', default=0, type=float, help='Angle of wedge formed between line'
                        'extending from pore center to monomer and line from pore center to vertically adjacent monomer'
                                                                                   'head group.')
    parser.add_argument('-angles', '--angles', nargs='+', default=[90, 90, 60], type=float, help='Angles between'
                        'box vectors')
    parser.add_argument('-seed', '--random_seed', default=False, type=int, help='Pass an integer to give a seed for '
                                                                                'random column displacement')

    # flags unique to equil.sh
    parser.add_argument('-ring_restraints', '--ring_restraints', default=["C", "C1", "C2", "C3", "C4", "C5"], nargs='+',
                        help='Name of atoms used to restrain head groups during initial equilibration.')
    parser.add_argument('-forces', default=[1000000, 3162, 56, 8, 3, 2, 1, 0], help='Sequence of force constants to'
                                                                                    'apply to ring restraints')
    parser.add_argument('-l_nvt', '--length_nvt', default=50, type=int, help='Length of restrained NVT simulations '
                                                                             '(ps)')
    parser.add_argument('-lb', '--length_berendsen', default=5000, type=int, help='Length of simulation using berendsen'
                                                                                  'pressure control')
    parser.add_argument('-lpr', '--length_Parrinello_Rahman', default=400000, type=int, help='Length of simulation '
                        'using Parrinello-Rahman pressure control')
    parser.add_argument('--restraint_residue', default='HII', type=str, help='Name of residue to which ring_restraint'
                                                                             'atoms belong')
    parser.add_argument('--restraint_axis', default='xyz', type=str, help='Axes along which to apply position '
                                                                          'restraints')

    return parser


def build(build_monomer, out, mon_per_col, ncol, radius, p2p, dbwl, pd, nopores=4, seed=False):

    # update once build.py is more class-based

    build = 'build.py -b %s -m %s -c %s -r %s -p %s -n %s -d %s -pd %s -o %s' % (build_monomer, mon_per_col, ncol,
                                                                                 radius, p2p, nopores, dbwl, pd, out)
    if seed:
        build += ' -seed %s' % seed

    subprocess.Popen(build.split()).wait()


def generate_input_files(gro, ensemble, length, barostat='berendsen', genvel=True, xlink=False, restraints=False,
                         frames=50, dt=0.002):

    # mostly uses defaults for now
    nstout = int(length / (dt * frames))
    mdp = SimulationMdp(gro, length=length, barostat=barostat, genvel=genvel, restraints=restraints, xlink=xlink,
                        nstxout=nstout, nstvout=nstout, nstfout=nstout, nstenergy=nstout)

    mdp.write_em_mdp()  # write energy minimization .mdp without asking

    if ensemble == 'npt':
        mdp.write_npt_mdp()
    elif ensemble == 'nvt':
        mdp.write_nvt_mdp()

    top = SystemTopology(gro, restraints=restraints, xlink=xlink)
    top.write_top()


def restrain(gro, build_monomer, force, axis, restraint_atoms):
    # update once restrain.py is class-based

    cmd = 'restrain.py -f %s -A %s -g %s -m %s -a %s' % (' '.join([str(force) for a in range(len(axis))]), axis, gro,
        build_monomer.split('.')[0], ' '.join(map(str, restraint_atoms)))

    p = subprocess.Popen(cmd.split())
    p.wait()


def simulate(mdp, top, config, out, restrained=False, mpi=False, np=4):

    if mpi:
        gmx = "mpirun -np %s gmx_mpi" % np
    else:
        gmx = "gmx"

    grompp = "%s grompp -f %s -p %s -c %s -o %s" % (gmx, mdp, top, config, out)

    if restrained:
        grompp += ' -r %s' % config

    p = subprocess.Popen(grompp.split())
    p.wait()

    mdrun = "%s mdrun -v -deffnm %s" % (gmx, out)
    p = subprocess.Popen(mdrun.split())
    p.wait()


def check_energy(logname='em.log'):

    nrg = 1

    with open('%s' % logname) as f:
        for line in f:
            if line.count("Potential Energy") == 1:
                nrg = float(line.split()[3])

    return nrg


if __name__ == "__main__":

    args = initialize().parse_args()

    os.environ["GMX_MAXBACKUP"] = "-1"  # stop GROMACS from making backups

    # initial build
    build(args.build_monomer, args.initial, args.monomers_per_column, args.ncolumns, args.pore_radius, args.p2p,
          args.dbwl, args.parallel_displaced, nopores=args.nopores, seed=args.random_seed)

    # generate input files once
    restrain(args.initial, args.build_monomer, args.forces[0], args.restraint_axis, args.ring_restraints)
    generate_input_files(args.initial, 'nvt', args.length_nvt, restraints=args.restraint_residue)

    simulate('em.mdp', 'topol.top', '%s' % args.initial, 'em', mpi=args.mpi, np=args.nproc, restrained=True)

    nrg = check_energy(logname='em.log')

    # rebuild until energy minimization reaches a negative potential energy
    while nrg > 0:

        build(args.build_monomer, args.initial, args.monomers_per_column, args.ncolumns, args.pore_radius, args.p2p,
                  args.dbwl, args.parallel_displaced, nopores=args.nopores)

        simulate('em.mdp', 'topol.top', '%s' % args.initial, 'em', mpi=args.mpi, np=args.nproc, restrained=True)

        nrg = check_energy(logname='em.log')

    simulate('nvt.mdp', 'topol.top', 'em.gro', 'nvt', mpi=args.mpi, np=args.nproc, restrained=True)

    copy = "cp nvt.gro %s.gro" % args.forces[0]
    p = subprocess.Popen(copy.split())
    copy = "cp nvt.trr %s.trr" % args.forces[0]
    p = subprocess.Popen(copy.split())

    generate_input_files('nvt.gro', 'nvt', args.length_nvt, genvel=False, restraints=args.restraint_residue)

    for f in args.forces[1:]:

        restrain(args.initial, args.build_monomer, f, args.restraint_axis, args.ring_restraints)
        simulate('nvt.mdp', 'topol.top', 'em.gro', 'nvt', mpi=args.mpi, np=args.nproc, restrained=True)

        copy = "cp nvt.gro %s.gro" % f
        p = subprocess.Popen(copy.split())
        p.wait()

        copy = "cp nvt.trr %s.trr" % f
        p = subprocess.Popen(copy.split())
        p.wait()

    generate_input_files(args.initial, 'npt', args.length_berendsen, genvel=False, barostat='berendsen')
    simulate('npt.mdp', 'topol.top', 'nvt.gro', 'berendsen', mpi=args.mpi, np=args.nproc)

    generate_input_files(args.initial, 'npt', args.length_Parrinello_Rahman, genvel=False, barostat='Parrinello-Rahman', frames=400)
    simulate('npt.mdp', 'topol.top', 'berendsen.gro', 'PR', mpi=args.mpi, np=args.nproc)
