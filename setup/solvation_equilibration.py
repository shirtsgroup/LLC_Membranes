#!/usr/bin/env python

import argparse
import subprocess
import sqlite3 as sql
from LLC_Membranes.llclib import topology, file_rw
from LLC_Membranes.analysis import solute_partitioning
from LLC_Membranes.setup import lc_class
import numpy as np
import os
import mdtraj as md

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in


def initialize():

    parser = argparse.ArgumentParser(description='Solvate system and adjust so there is a certain wt % of water')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file to add water to')
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

    # flags unique to equil.sh
    parser.add_argument('-ring_restraints', '--ring_restraints', default=["C", "C1", "C2", "C3", "C4", "C5"], nargs='+',
                        help='Name of atoms used to restrain head groups during initial equilibration.')

    return parser


class System(object):

    def __init__(self, args, solute='HOH'):
        """

        :param args: arguments fed to argparse
        :param solute: name of solute used to solvate system
        """

        self.args = args

        # get build monomer molecular weight and calculate mw of entire dry system
        self.build_monomer = topology.Residue(args.build_monomer.split('.')[0])
        if self.build_monomer.residues:
            self.build_monomer_mw = sum([x.mw for x in self.build_monomer.residues])
        else:
            self.build_monomer_mw = self.build_monomer.mw

        nmon = args.nopores * args.ncolumns * args.monomers_per_column  # number of build monomers in the system
        self.dry_mass = nmon * self.build_monomer_mw  # mass of dry system

        # calculate required water in pores and tails
        self.water = topology.Residue(solute)
        self.args.weight_percent /= 100.0  # convert to fraction
        self.total_water = int((self.args.weight_percent * nmon * self.build_monomer_mw) / (self.water.mw *
                            (1 - self.args.weight_percent)))
        tail_water = self.total_water / (self.args.ratio + 1)
        self.pore_water = int(args.ratio * tail_water)
        self.tail_water = int(tail_water)

        self.r = 0  # pore radius
        self.converged = False
        self.solvated = None  # an object that will describe solvated systems

    def query_database(self, database='water_content.db'):

        # read database of pore radii and associated water contents
        connection = sql.connect("%s/%s" % (location, database))  # database created in this directory
        crsr = connection.cursor()
        sql_command = "select nwater, radius from radii where monomer = '%s' and pd_angle = %.2f and dbwl = %.2f and " \
                      "mon_per_col = %d;" % (self.args.build_monomer.split('.')[0], self.args.parallel_displaced,
                                             self.args.dbwl, self.args.monomers_per_column)

        try:
            sql_output = crsr.execute(sql_command).fetchall()
        except sql.OperationalError:
            sql_output = []

        if sql_output:

            # assumes nwater scales with radii (which might be erroneous for data points that are close together)
            nwater = sorted([x[0] for x in sql_output])
            radii = sorted([float(x[1]) for x in sql_output])
            print(radii, nwater, self.pore_water)

            if len(nwater) == 1:

                if self.pore_water > nwater[0]:
                    self.r = radii[0] + self.args.guess_stride
                elif self.pore_water < nwater[0]:
                    self.r = radii[0] - self.args.guess_stride
                else:
                    self.r = radii[0]

            elif self.pore_water < min(nwater):

                self.r = min(radii) - self.args.guess_stride

            elif self.pore_water > max(nwater):

                self.r = max(radii) + self.args.guess_stride

            else:

                bin = np.digitize(self.pore_water, nwater)

                upper_bound, lower_bound = nwater[bin], nwater[
                    bin - 1]  # upper bound is exclusive. Lower bound is inclusive

                # calculate water content if an entry doesn't already exist
                if abs(self.pore_water - lower_bound) > self.args.tolerance:
                    # linearly interpolate what the next pore radius should be based on desired amount of water in pore
                    interpolation = (self.pore_water - lower_bound) / (upper_bound - lower_bound)
                    self.r = radii[bin - 1] + (radii[bin] - radii[bin - 1]) * interpolation
                else:
                    self.r = radii[bin - 1]
        else:

            self.r = (self.args.guess_range[1] - self.args.guess_range[0]) * np.random.sample() + self.args.guess_range[0]

        connection.close()

    def equilibrate(self):

        subprocess.call(['build.py', '-b', '%s' % self.args.build_monomer, '-m', '%s' % self.args.monomers_per_column,
                         '-c', '%s' % self.args.ncolumns, '-r', '%s' % self.r, '-p', '%s' % self.args.p2p, '-n', '%s'
                         % self.args.nopores, '-d', '%s' % self.args.dbwl, '-pd', '%s' % self.args.parallel_displaced])

        subprocess.call(['equil.sh', '-q', '1', '-b', '%s' % self.args.build_monomer.split('.')[0], '-r',
                         '%s' % ' '.join(self.args.ring_restraints), '-m', '%s' % self.args.mpi, '-p', '%s'
                         % self.args.nproc])

    def calculate_pore_water(self):

        # solvate the system
        if self.args.mpi:
            gmx = "mpirun -np %s gmx_mpi" % self.args.np
        else:
            gmx = "gmx"

        cmd = "%s solvate -cp 1000000.gro -cs spc216.gro -o solvated.gro -p topol.top" % gmx
        subprocess.call(cmd.split())

        pore_defining_atoms = lc_class.LC(self.args.build_monomer).pore_defining_atoms
        self.solvated = solute_partitioning.System('solvated.gro', pore_defining_atoms, 'SOL')
        self.solvated.locate_pore_centers()

        # radius based on reference atom of lc_class, but partition based on pore_defining_atoms. Need to make choice or leave it
        self.solvated.partition(self.r)

        if abs(self.pore_water - len(self.solvated.pore_water[0])) <= self.args.tolerance:
            self.converged = True

        if self.pore_water != self.solvated.pore_water[0]:  # avoid duplicates
            self.update_database()

    def update_database(self, database='water_content.db'):

        connection = sql.connect("%s/%s" % (location, database))  # database created in this directory

        crsr = connection.cursor()

        sql_command = "insert into radii (monomer, radius, mon_per_col, nwater, pd_angle, dbwl) values ('%s', %.3f, %d," \
                      "%d, %.2f, %.2f);" % (self.args.build_monomer.split('.')[0], self.r, self.args.monomers_per_column,
                        len(self.solvated.pore_water[0]), self.args.parallel_displaced, self.args.dbwl)

        crsr.execute(sql_command)

        connection.commit()
        connection.close()

    def write_final_configuration(self):

        # only do this if we have converged on the correct number of waters
        water_indices = []
        for i in self.solvated.tail_water[0]:  # tail_water indices are given as the center of mass of each water
            water_indices += self.solvated.residue_indices[(self.water.natoms * i): self.water.natoms * (i + 1)].tolist()

        keep = np.full(self.solvated.pos[1], True, dtype=bool)  # array of True. Booleans are fast
        keep[water_indices] = False  # set pore indices to False

        file_rw.write_gro_pos(self.solvated.pos[0, keep, :], self.args.output, ids=np.array(self.solvated.ids)[keep],
                              res=np.array(self.solvated.res)[keep])


if __name__ == "__main__":

    args = initialize().parse_args()

    sys = System(args)

    while not sys.converged:
        sys.query_database()
        print(sys.r)
        sys.equilibrate()
        sys.calculate_pore_water()

    sys.write_final_configuration()