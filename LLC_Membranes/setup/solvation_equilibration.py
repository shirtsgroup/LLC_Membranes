#!/usr/bin/env python

import argparse
import subprocess
import sqlite3 as sql
from LLC_Membranes.llclib import topology, file_rw
from LLC_Membranes.analysis import solute_partitioning
from LLC_Membranes.setup import lc_class, equil, solvate_tails
import mdtraj as md
import numpy as np
import os


def initialize():

    parser = argparse.ArgumentParser(description='Solvate system and adjust so there is a certain wt % of water')

    parser.add_argument('-ratio', '--ratio', default=2, type=float, help='Ratio of water in pores to water in tails')
    parser.add_argument('-wt', '--weight_percent', default=10, type=float, help='Total weight percent of water')
    parser.add_argument('-tol', '--tolerance', default=1, type=int, help='Number of water molecules')
    parser.add_argument('-guess_range', default=[.4, 1], nargs='+', help='If water_content.db has no entries for the '
                        'build monomer, an initial radius will be randomly selected from this range')
    parser.add_argument('-guess_stride', default=0.2, type=float, help='How far above/below the highest/lowest value to'
                        'make the next guess at pore radius if you need more/less water than the bounds of the water '
                        'content database (nm)')
    parser.add_argument('-o', '--output', default='solvated_final.gro', help='Name of fully solvated output file')
    parser.add_argument('-seed', '--random_seed', default=0, type=int, help='Numpy random seed')

    # parallelization
    parser.add_argument('-mpi', '--mpi', action="store_true", help='Run MD simulations in parallel')
    parser.add_argument('-np', '--nproc', default=4, help='Number of MPI processes')

    # same flags as to build.py
    parser.add_argument('-b', '--build_monomer', default='NAcarb11V', type=str, help='Name of single monomer used to '
                                                                                     'build full system')
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

    # flags unique to equil.py
    parser.add_argument('-ring_restraints', '--ring_restraints', default=["C", "C1", "C2", "C3", "C4", "C5"], nargs='+',
                        help='Name of atoms used to restrain head groups during initial equilibration.')
    parser.add_argument('-forces', default=[1000000, 3162, 56, 8, 3, 2, 1, 0], help='Sequence of force constants to'
                                                                                    'apply to ring restraints')
    parser.add_argument('--restraint_residue', default='HII', type=str, help='Name of residue to which ring_restraint'
                                                                             'atoms belong')
    parser.add_argument('--restraint_axis', default='xyz', type=str, help='Axes along which to apply position '
                                                                          'restraints')
    parser.add_argument('-l_nvt', '--length_nvt', default=50, type=int, help='Length of restrained NVT simulations '
                                                                             '(ps)')
    parser.add_argument('-lb', '--length_berendsen', default=5000, type=int, help='Length of simulation using berendsen'
                                                                                  'pressure control')
    parser.add_argument('-lpr', '--length_Parrinello_Rahman', default=400000, type=int, help='Length of simulation '
                        'using Parrinello-Rahman pressure control')

    return parser


class System(object):

    def __init__(self, build_monomer, weight_percent, ratio, tolerance=1, solute='HOH', nopores=4, ncolumns=5,
                 monomers_per_column=20, p2p=4.5, parallel_displaced=0, dbwl=0.37, random_seed=0, mpi=False, nproc=4):
        """ Get unit cell build parameters and determine the number of water molecules required for each region in order
        to have the correct final composition.

        :param build_monomer: Name of liquid crystal monomer to use to build unit cell (no file extension)
        :param weight_percent: Weight percent of water in the whole system
        :param ratio: Ratio of water in the pores to water in the tails
        :param tolerance: Acceptable error in the total number of water molecules placed in the pore region
        :param solute: name of solute used to solvate system
        :param nopores: Number of pores in unit cell
        :param ncolumns: Number of stacked monomer columns per pore
        :param monomers_per_column: Number of monomers stacked into a column
        :param p2p: Distance between pores (nm)
        :param parallel_displaced: Angle of wedge formed between line extending from pore center to monomer and line \
        from pore center to vertically adjacent monomer head group
        :param dbwl: Distance between stacked monomers (nm)
        :param random_seed: Monomer columns are randomly displaced in the z-direction when the initial configuration \
        is built. Specify a random seed so that the same displacement is achieved each time the radius is changed. I \
        think this helps with convergence, but it hasn't been extensively tested.
        :param mpi: Run the MPI version of GROMACS
        :param np: number of MPI process if mpi = True

        :type build_monomer: str
        :type weight_percent: float
        :type ratio: float
        :type tolerance: int
        :type solute: str
        :type nopores: int
        :type ncolumns: int
        :type monomers_per_column: int
        :type p2p: float
        :type parallel_displaced: float
        :type dbwl: float
        :type random_seed: int
        :type mpi: bool
        :type np: int
        """

        # Initialize variables needed later
        self.location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        self.parallel_displaced = parallel_displaced
        self.dbwl = dbwl
        self.random_seed = random_seed
        self.monomers_per_column = monomers_per_column
        self.ncolumns = ncolumns
        self.p2p = p2p
        self.nopores = nopores
        self.mpi = mpi
        self.np = nproc
        self.tolerance = tolerance
        self.build_monomer = topology.LC(build_monomer)

        t = md.load('%s/../top/topologies/%s.gro' % (self.location, build_monomer))
        residues = set([a.residue.name for a in t.topology.atoms])

        # since liquid crystals with counterions contain multiple residues
        self.build_monomer_mw = 0
        for r in residues:
            self.build_monomer_mw += topology.Residue(r).MW

        nmon = self.nopores * self.ncolumns * self.monomers_per_column  # number of build monomers in the system
        self.dry_mass = nmon * self.build_monomer_mw  # mass of dry system

        # calculate required water in pores and tails
        self.water = topology.Residue(solute)
        self.weight_percent = weight_percent / 100.0  # convert to fraction
        self.total_water = int((self.weight_percent * nmon * self.build_monomer_mw) / (self.water.MW *
                            (1 - self.weight_percent)))
        tail_water = self.total_water / (ratio + 1)
        self.pore_water = int(ratio * tail_water)
        self.tail_water = int(tail_water)

        self.r = 0  # pore radius
        self.converged = False
        self.solvated = None  # an object that will describe solvated systems

    def query_database(self, database='water_content.db', guess_range=(0.4, 1), guess_stride=0.2):
        """ Read an SQL database, of number of water molecules versus pore radius, to inform the next choice of pore
        radius. The decision process bisects the two radii with water contents closest to the desired concentration. If
        there is a single data point, an educated guess is made. If there are no data points, then make a random guess.

        :param database: Name of database
        :param guess_range: If water_content.db has no entries for the build monomer, an initial radius will be \
        randomly selected from this range
        :param guess_stride: How far above/below the highest/lowest value to make the next guess at pore radius if you \
        need more/less water than the bounds of the water

        :type database: str
        :type guess_range: tuple
        :type guess_stride: float
        """

        # read database of pore radii and associated water contents
        connection = sql.connect("%s/%s" % (self.location, database))  # database created in this directory
        crsr = connection.cursor()
        sql_command = "select nwater, radius from radii where monomer = '%s' and pd_angle = %.2f and dbwl = %.2f and " \
                      "mon_per_col = %d and seed = %s;" % (self.build_monomer.name, self.parallel_displaced, self.dbwl,
                                                           self.monomers_per_column, self.random_seed)

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
                    self.r = radii[0] + guess_stride
                elif self.pore_water < nwater[0]:
                    self.r = radii[0] - guess_stride
                else:
                    self.r = radii[0]

            elif self.pore_water < min(nwater):

                self.r = min(radii) - guess_stride

            elif self.pore_water > max(nwater):

                self.r = max(radii) + guess_stride

            else:

                bin = np.digitize(self.pore_water, nwater)

                upper_bound, lower_bound = nwater[bin], nwater[
                    bin - 1]  # upper bound is exclusive. Lower bound is inclusive

                # calculate water content if an entry doesn't already exist
                if abs(self.pore_water - lower_bound) > self.tolerance:
                    # linearly interpolate what the next pore radius should be based on desired amount of water in pore
                    interpolation = (self.pore_water - lower_bound) / (upper_bound - lower_bound)
                    self.r = radii[bin - 1] + (radii[bin] - radii[bin - 1]) * interpolation
                else:
                    self.r = radii[bin - 1]
        else:

            self.r = (guess_range[1] - guess_range[0]) * np.random.sample() + guess_range[0]

        connection.close()

    def build(self):

        equil.build('%s.gro' % self.build_monomer.name, 'initial.gro', self.monomers_per_column, self.ncolumns,
                    self.r, self.p2p, self.dbwl, self.parallel_displaced,
                    nopores=self.nopores, seed=self.random_seed)

    def restrain(self, force, restraint_axis='xyz', ring_restraints=("C", "C1", "C2", "C3", "C4", "C5")):

        equil.restrain('initial.gro', self.build_monomer.name, force, restraint_axis, ring_restraints)

    def input_files(self, gro, ensemble, length=50, restraint_residue='HII'):

        equil.generate_input_files(gro, ensemble, length, restraints=restraint_residue)

    def put_in_box(self, gro, tpr):

        p1 = subprocess.Popen(['echo', '0'], stdout=subprocess.PIPE)
        trjconv = "gmx trjconv -f %s -o %s -pbc atom -s %s -ur tric" % (gro, gro, tpr)
        p2 = subprocess.Popen(trjconv.split(), stdin=p1.stdout)
        p2.wait()

    def equilibrate(self, input_files=True, length=50, force=1000000, restraint_residue='HII', restraint_axis='xyz',
                    ring_restraints=("C", "C1", "C2", "C3", "C4", "C5")):
        """ Simulate the unit cell with restraints placed on the head group

        :param input_files: Generate GROMACS .mdp and topology files
        :param length: Simulation length (ps)
        :param force: Force with which to restrain head groups kJ mol^-1 nm^-2
        :param restraint_residue: Name of residue to which position restraints will be applied
        :param restraint_axis: Axis/axes along which head groups should be restrained
        :param ring_restraints: Names of head group atoms

        :type input_files: bool
        :type length: int
        :type force: float or int
        :type restraint_residue: str
        :type restraint_axis: str
        :type ring_restraints: tuple
        """

        self.build()  # build initial configuration

        if input_files:  # doesn't need to be done every time
            self.restrain(force, restraint_axis=restraint_axis, ring_restraints=ring_restraints)
            self.input_files('initial.gro', 'nvt', length=length, restraint_residue=restraint_residue)

        nrg = 1
        while nrg > 0:

            self.build()

            equil.simulate('em.mdp', 'topol.top', 'initial.gro', 'em', mpi=self.mpi, np=self.np, restrained=True)

            nrg = equil.check_energy(logname='em.log')

            if nrg > 0:  # choose a new random seed if energy minimization doesn't work
                self.random_seed = np.random.randint(0, 4294967295)

        cp = 'cp em.gro %s.gro' % force
        p = subprocess.Popen(cp.split())
        p.wait()

        self.put_in_box('%s.gro' % force, 'em.tpr')

    def calculate_pore_water(self, config):
        """ Determine the total number of water molecules within the pore radius. Update database with system
        configuration parameters.

        :param config: Name of .gro configuration of which to calculate pore water content

        :type config: str
        """

        # solvate the system
        if self.mpi:
            gmx = "mpirun -np %s gmx_mpi" % self.np
        else:
            gmx = "gmx"

        cmd = "%s solvate -cp %s.gro -cs spc216.gro -o solvated.gro -p topol.top" % (gmx, config)
        subprocess.call(cmd.split())

        self.put_in_box('solvated.gro', 'solvated.gro')

        self.solvated = solute_partitioning.System('solvated.gro', self.build_monomer.name, 'SOL')
        self.solvated.locate_pore_centers()

        # radius based on reference atom, but partition based on pore_defining_atoms. Need to make choice or leave it
        self.solvated.partition(self.r)

        if abs(self.pore_water - len(self.solvated.pore_water[0])) <= self.tolerance:
            self.converged = True

        if self.pore_water != self.solvated.pore_water[0]:  # avoid duplicates
            self.update_database()

    def update_database(self, database='water_content.db'):

        connection = sql.connect("%s/%s" % (self.location, database))  # database created in this directory

        crsr = connection.cursor()

        sql_command = "insert into radii (monomer, radius, mon_per_col, nwater, pd_angle, dbwl, seed) values ('%s'," \
                      "%.6f, %d, %d, %.2f, %.2f, %d);" % (self.build_monomer.name, self.r, self.monomers_per_column,
                                                          len(self.solvated.pore_water[0]), self.parallel_displaced,
                                                          self.dbwl, self.random_seed)

        crsr.execute(sql_command)

        connection.commit()
        connection.close()

    def write_final_pore_configuration(self):
        """ Write the configuration with the correct number of water molecules placed in the pore region, then
        rwrite topology files.
        """

        # only do this if we have converged on the correct number of waters
        water_indices = []

        for i in self.solvated.tail_water[0]:  # tail_water indices are given as the center of mass of each water
            water_indices += self.solvated.residue_indices[(self.water.natoms * i): self.water.natoms * (i + 1)].tolist()

        keep = np.full(self.solvated.pos.shape[1], True, dtype=bool)  # array of True. Booleans are fast
        keep[water_indices] = False  # set pore indices to False

        # change all 'HOH' to 'SOL' because mdtraj changed it
        res = np.array(self.solvated.res)
        res[np.where(np.array(self.solvated.res) == 'HOH')[0]] = 'SOL'

        file_rw.write_gro_pos(self.solvated.pos[0, keep, :], 'solvated_pores.gro', ids=np.array(self.solvated.ids)[keep],
                              res=res[keep], box=self.solvated.box)

        # rewrite topology files
        self.input_files('solvated_pores.gro', 'nvt')

    def place_water_tails(self, output):
        """ Place the desired number of water molecules in the tails. This is done by randomly inserting water molecules
        in close proximity to the tails, one-by-one. A short energy minimzation is performed between each insertion.

        :param output: Name of final configuration

        :type output: str
        """

        tails = solvate_tails.System('solvated_pores.gro', 'topol.top', self.build_monomer.name, rbounds=[0.3, 1],
                                     restraints=True, mpi=self.mpi, nproc=self.np)

        tails.insert_all_water(self.tail_water, output=output, final_topname='topol.top')

    def full_equilibration(self, forces, fully_solvated='solvated_final.gro', l_nvt=50, l_berendsen=5000, l_pr=400000,
                           restraint_residue='HII', restraint_axis='xyz',
                           ring_restraints=("C", "C1", "C2", "C3", "C4", "C5")):
        """ Simulate the unit cell with a sequence of decreasing restraints placed on the head group

        :param forces: sequence of forces to apply to ring_restraints (:math:`\dfrac{kJ}{mol~nm^2}`)
        :param fully_solvated: name of fully solvated coordinate file
        :param l_nvt: Length of short restrained simulations (ps)
        :param l_berendsen: Length of equilibration simulation run with berendsen pressure control (ps)
        :param l_pr: Length of long equilibration simulation run with Parrinello-Rahman pressure control (ps)
        :param restraint_residue: Name of residue to which position restraints will be applied
        :param restraint_axis: Axis/axes along which head groups should be restrained
        :param ring_restraints: Names of head group atoms

        :type forces: tuple or list
        :type fully_solvated: str
        :type l_nvt: int
        :type l_berendsen: int
        :type l_pr: int
        :type restraint_residue: str
        :type restraint_axis: str
        :type ring_restraints: tuple
        """

        equil.simulate('nvt.mdp', 'topol.top', '%s' % fully_solvated, 'nvt', mpi=self.mpi, np=self.np, restrained=True)

        cp = "cp nvt.gro em.gro"
        subprocess.Popen(cp.split()).wait()

        cp = "cp nvt.gro %s.gro" % forces[0]
        subprocess.Popen(cp.split()).wait()

        equil.generate_input_files('nvt.gro', 'nvt', l_nvt, genvel=False, restraints=restraint_residue)

        for f in forces[1:]:

            equil.restrain(fully_solvated, self.build_monomer.name, f, restraint_axis, ring_restraints)
            equil.simulate('nvt.mdp', 'topol.top', 'em.gro', 'nvt', mpi=self.mpi, np=self.np, restrained=True)

            cp = "cp nvt.gro %s.gro" % f
            subprocess.Popen(cp.split()).wait()

            cp = "cp nvt.trr %s.trr" % f
            subprocess.Popen(cp.split()).wait()

        equil.generate_input_files(fully_solvated, 'npt', l_berendsen, genvel=False, barostat='berendsen', frames=50)
        equil.simulate('npt.mdp', 'topol.top', 'nvt.gro', 'berendsen', mpi=self.mpi, np=self.np)

        equil.generate_input_files('berendsen.gro', 'npt', l_pr, genvel=False, barostat='Parrinello-Rahman', frames=500)
        equil.simulate('npt.mdp', 'topol.top', 'berendsen.gro', 'PR', mpi=self.mpi, np=self.np)


if __name__ == "__main__":

    os.environ["GMX_MAXBACKUP"] = "-1"  # stop GROMACS from making backups

    args = initialize().parse_args()

    sys = System(args.build_monomer, args.weight_percent, args.ratio, solute='HOH', nopores=args.nopores,
                 ncolumns=args.ncolumns, monomers_per_column=args.monomers_per_column, p2p=args.p2p,
                 parallel_displaced=args.parallel_displaced, dbwl=args.dbwl, random_seed=args.random_seed,
                 mpi=args.mpi, nproc=args.nproc, tolerance=args.tolerance)

    while not sys.converged:
        sys.query_database(database='water_content.db', guess_range=args.guess_range, guess_stride=args.guess_stride)
        sys.equilibrate(input_files=True, length=args.length_nvt, force=args.forces[0],
                        restraint_residue=args.restraint_residue, restraint_axis=args.restraint_axis,
                        ring_restraints=args.ring_restraints)
        sys.calculate_pore_water(args.forces[0])

    sys.write_final_pore_configuration()
    sys.place_water_tails(args.output)
    sys.full_equilibration(args.forces, fully_solvated=args.output, l_berendsen=args.length_berendsen,
                           l_nvt=args.length_nvt, l_pr=args.length_Parrinello_Rahman,
                           restraint_axis=args.restraint_axis, restraint_residue=args.restraint_residue,
                           ring_restraints=args.ring_restraints)
