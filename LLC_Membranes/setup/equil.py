#!/usr/bin/env python

import argparse
import os
import subprocess
from LLC_Membranes.setup.gentop import SystemTopology
from LLC_Membranes.setup.genmdp import SimulationMdp
from LLC_Membranes.setup.restrain import RestrainedTopology
from LLC_Membranes.llclib import topology, physical, file_rw, gromacs
import mdtraj as md
import numpy as np

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
    parser.add_argument('-b', '--build_monomer', default='NAcarb11V.gro', nargs='+', type=str, help='Name of single '
                        'monomer structure file (.gro format) used to build full system')
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
    parser.add_argument('-mf', '--mol_frac', nargs='+', default=1., type=float, help='If using the -random flag, this gives'
                        'the relative amount of each monomer to add. List the fractions in the same order as'
                        '-build_monomer')

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
    parser.add_argument('--restraint_residue', default='HII', nargs='+', type=str,
                        help='Name of residue to which ring_restraint atoms belong')
    parser.add_argument('--restraint_axis', default='xyz', type=str, help='Axes along which to apply position '
                                                                          'restraints')

    parser.add_argument('--continue_initial_config', default=None, help='Equilibrate energy minimized initial '
                                                                        'configuration.')
    parser.add_argument('-try', '--ntries', default=10, type=int, help='The number of times to try to energy minimize an initial'
                        'configuration before resorting to scaling + shrinking it.')

    return parser

# TODO: make a yaml


class HexagonalPhaseEquilibration:

    def __init__(self, build_monomers, mpi=False, nprocesses=4):

        self.lc = [topology.LC(m) for m in build_monomers]

        self.gro_name = None
        self.top = None
        self.mdp = None
        self.mpi = mpi
        self.nprocesses = nprocesses

    def build(self, build_monomer, out, mon_per_col, ncol, radius, p2p, dbwl, pd, nopores=4, seed=False,
              mole_fraction=1.):
        # TODO: incorporate class-based functionality

        # update once build.py is more class-based

        build = 'build.py -phase h2 -m %s -c %s -r %s -p %s -n %s -d %s -pd %s -o %s -b' % (mon_per_col, ncol,
                radius, p2p, nopores, dbwl, pd, out)

        for mon in build_monomer:
            build += ' %s' % mon

        if type(mole_fraction) is float:
            mole_fraction = [mole_fraction]

        build += ' -mf'
        for mf in mole_fraction:
            build += ' %s' % mf

        if seed:
            build += ' -seed %s' % seed

        subprocess.Popen(build.split()).wait()

        self.gro_name = out

    def generate_input_files(self, ensemble, length, barostat='berendsen', genvel=True, xlink=False, restraints=False,
                             frames=50, dt=0.002):

        # mostly uses defaults for now
        nstout = int(length / (dt * frames))
        self.mdp = SimulationMdp(self.gro_name, length=length, barostat=barostat, genvel=genvel, restraints=restraints,
                                 xlink=xlink, nstxout=nstout, nstvout=nstout, nstfout=nstout, nstenergy=nstout)

        self.mdp.write_em_mdp()  # write energy minimization .mdp without asking

        if ensemble == 'npt':
            self.mdp.write_npt_mdp()
        elif ensemble == 'nvt':
            self.mdp.write_nvt_mdp()

        self.top = SystemTopology(self.gro_name, restraints=restraints, xlink=xlink)
        self.top.write_top()

    def restrain(self, build_monomer, force, axis, restraint_atoms):

        top = RestrainedTopology(self.gro_name, build_monomer, restraint_atoms, com=False)
        top.add_position_restraints(axis, [force for _ in range(len(axis))])
        top.write_topology()

    def shrink_unit_cell(self, start, stop=1, step=0.2):
        """ Sequentially energy minimize a structure to a packing structure by isotropically scaling the distance between
        monomer head groups

        :return:
        """

        for i in np.linspace(start - step, stop, int(round((start - stop) / step))):

            f = i / (i + step)
            self.scale_columns(f)
            print('Attempting energy minimization with box lengths %.4f times desired size' % i, end='', flush=True)
            nrg = gromacs.simulate(self.mdp.em_mdp_name, self.top.name, self.gro_name,
                                   self.mdp.em_mdp_name.split('.')[0], em_energy=True, verbose=False, mpi=self.mpi,
                                   nprocesses=self.nprocesses, restraints=True)

            if nrg >= 0:
                print('...Failed. Will shrink slower on next iteration.')
                self.gro_name = 'scaled_%.4f.gro' % (i + step)
                self.shrink_unit_cell(i + step, i, step / 2)  # reduce shrink rate
            else:
                print('Running short NVT simulation...', end='', flush=True)
                self.gro_name = self.mdp.em_mdp_name.split('.')[0] + '.gro'
                self.simulate(self.mdp.nvt_mdp_name, self.top.name, 'nvt_%.4f' % i, restrained=True)
                print('...Success!')

            cp = 'cp em.gro scaled_%.4f.gro' % i
            p = subprocess.Popen(cp.split())
            p.wait()

    def scale_columns(self, scale_factor):

        pore_defining_atoms = [self.lc[i].pore_defining_atoms for i in range(len(self.lc))]

        t = md.load(self.gro_name)
        resnames = [a.residue.name for a in t.topology.atoms]
        ids = [a.name for a in t.topology.atoms]

        indices = [[] for _ in range(len(self.lc))]
        monomers = {l.name: i for i, l in enumerate(self.lc)}

        for a in t.topology.atoms:
            try:
                mon = monomers[a.residue.name]
                if a.name in pore_defining_atoms[mon]:
                    indices[mon].append(a.index)
            except KeyError:
                pass

        mass = [[self.lc[i].mass[a] for a in self.lc[i].LC_names if a in pda]
                for i, pda in enumerate(pore_defining_atoms)]

        com = []
        for m in range(len(monomers.values())):
            com.append(physical.center_of_mass(t.xyz[:, indices[m], :], mass[m]))

        # rebuild unit cell with scaled center-of-mass z coordinates
        nres = [len(indices[i]) // len(pore_defining_atoms[i]) for i in range(len(self.lc))]
        ndx = 0
        n = [0, 0]
        for i in range(sum(nres)):
            mon = monomers[resnames[ndx]]
            end = ndx + self.lc[mon].natoms
            ref_z = com[mon][0, n[mon], 2]
            scaled_com = ref_z * scale_factor
            t.xyz[0, ndx:end, 2] += (scaled_com - ref_z)
            ndx = end  # - self.lc[mon].no_ions  -- haven't tested
            n[mon] += 1

        ucell = t.unitcell_vectors[0, ...]
        ucell[2, 2] *= scale_factor

        self.gro_name = 'scaled.gro'
        file_rw.write_gro_pos(t.xyz[0, ...], self.gro_name, ids=ids, res=resnames, ucell=ucell)

    def simulate(self, mdp, top, out, restrained=False):

        if self.mpi:
            gmx = "mpirun -np %s gmx_mpi" % self.nprocesses
        else:
            gmx = "gmx"

        grompp = "%s grompp -f %s -p %s -c %s -o %s" % (gmx, mdp, top, self.gro_name, out)

        if restrained:
            grompp += ' -r %s' % self.gro_name

        p = subprocess.Popen(grompp.split())
        p.wait()

        mdrun = "%s mdrun -v -deffnm %s" % (gmx, out)
        p = subprocess.Popen(mdrun.split())
        p.wait()

        self.gro_name = out + '.gro'


def check_energy(logname='em.log'):

    nrg = 1

    with open('%s' % logname) as f:
        for line in f:
            if line.count("Potential Energy") == 1:
                nrg = float(line.split()[3])

    return nrg


if __name__ == "__main__":

    # TODO: replace equil.simulate() with gromacs.simulate()
    args = initialize().parse_args()

    os.environ["GMX_MAXBACKUP"] = "-1"  # stop GROMACS from making backups

    equil = HexagonalPhaseEquilibration(args.build_monomer, mpi=args.mpi, nprocesses=args.nproc)

    lc = [topology.LC(mon) for mon in args.build_monomer]
    atoms = [l.pore_defining_atoms for l in lc]

    if args.continue_initial_config is None:

        # initial build
        equil.build(args.build_monomer, args.initial, args.monomers_per_column, args.ncolumns, args.pore_radius, args.p2p,
              args.dbwl, args.parallel_displaced, nopores=args.nopores, seed=args.random_seed, mole_fraction=args.mol_frac)

        # generate input files once
        equil.restrain(args.build_monomer, args.forces[0], args.restraint_axis, atoms)

        equil.generate_input_files('nvt', args.length_nvt, restraints=args.restraint_residue)

        # try energy minimizing
        nrg = gromacs.simulate(equil.mdp.em_mdp_name, equil.top.name, equil.gro_name, 'em', em_energy=True,
                               verbose=True, mpi=equil.mpi, nprocesses=4, restraints=True)

        try_ = 0  # try rebuilding the system ntries times to see if we get any working configs
        while nrg > 0 and try_ < args.ntries:  # rebuild if that doesn't work
            print('Try # %d' % (try_ + 1))
            equil.build(args.build_monomer, args.initial, args.monomers_per_column, args.ncolumns, args.pore_radius,
                        args.p2p,
                        args.dbwl, args.parallel_displaced, nopores=args.nopores, seed=args.random_seed,
                        mole_fraction=args.mol_frac)

            equil.restrain(args.build_monomer, args.forces[0], args.restraint_axis, atoms)

            nrg = gromacs.simulate('em.mdp', equil.top.name, equil.gro_name, 'em', em_energy=True,
                                   verbose=False, mpi=equil.mpi, nprocesses=4, restraints=True)

            try_ += 1

        if nrg > 0:  # if energy minimization still doesn't work

            equil.scale_columns(2)  # scale the unit cell

            # try energy minimization again
            nrg = gromacs.simulate('em.mdp', equil.top.name, equil.gro_name, 'em', em_energy=True,
                                   verbose=True, mpi=equil.mpi, nprocesses=4, restraints=True)

            while nrg > 0:  # rebuild if that doesn't work

                equil.build(args.build_monomer, args.initial, args.monomers_per_column, args.ncolumns, args.pore_radius,
                            args.p2p,
                            args.dbwl, args.parallel_displaced, nopores=args.nopores, seed=args.random_seed,
                            mole_fraction=args.mol_frac)

                if len(args.build_monomer) > 1:  # order of monomers changes if there is a mix. Otherwise can keep top
                    equil.restrain(args.build_monomer, args.forces[0], args.restraint_axis, atoms)

                equil.scale_columns(2)

                nrg = gromacs.simulate('em.mdp', equil.top.name, equil.gro_name, 'em', em_energy=True,
                                       verbose=True, mpi=equil.mpi, nprocesses=4, restraints=True)

            cp = 'cp em.gro scaled_2.0000.gro'
            p = subprocess.Popen(cp.split())
            p.wait()

            equil.shrink_unit_cell(2, 1, 0.2)

        equil.gro_name = 'em.gro'

    else:

        equil.gro_name = args.continue_initial_config

        equil.restrain(args.build_monomer, args.forces[0], args.restraint_axis, atoms)

        equil.generate_input_files('nvt', args.length_nvt, restraints=args.restraint_residue)

    equil.simulate('nvt.mdp', 'topol.top', 'nvt', restrained=True)

    copy = "cp nvt.gro %s.gro" % args.forces[0]
    p = subprocess.Popen(copy.split())
    copy = "cp nvt.trr %s.trr" % args.forces[0]
    p = subprocess.Popen(copy.split())

    equil.generate_input_files('nvt.gro', args.length_nvt, genvel=False, restraints=args.restraint_residue)

    exit()

    for f in args.forces[1:]:

        equil.restrain(args.build_monomer, f, args.restraint_axis, atoms)
        equil.simulate(equil.mdp.nvt_mdp_name, equil.top.name, equil.mdp.nvt_mdp_name.split('.')[0], restrained=True)

        copy = "cp nvt.gro %s.gro" % f
        p = subprocess.Popen(copy.split())
        p.wait()

        copy = "cp nvt.trr %s.trr" % f
        p = subprocess.Popen(copy.split())
        p.wait()

    equil.generate_input_files('npt', args.length_berendsen, genvel=False, barostat='berendsen')
    equil.simulate('npt.mdp', 'topol.top', 'berendsen')

    equil.generate_input_files('npt', args.length_Parrinello_Rahman, genvel=False, barostat='Parrinello-Rahman',
                               frames=400)
    equil.simulate('npt.mdp', 'topol.top', 'PR')
