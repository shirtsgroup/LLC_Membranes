#!/usr/bin/env python

import argparse
import place_solutes_pores
import numpy as np
import mdtraj as md
import subprocess
from gentop import SystemTopology
from genmdp import SimulationMdp


def initialize():
    parser = argparse.ArgumentParser(description='Set up and run umbrella sampling simulations')

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Equilibrated coordinate file')
    parser.add_argument('-s', '--solute', default='ETH', type=str, help='Name of solute residue to pull through pore')
    parser.add_argument('-equil_length', default=10000, type=int, help='Equilibration simulation length (ps)')
    parser.add_argument('-umbrella_length', default=10000, type=int, help='Length of each umbrella simulation')

    args = parser.parse_args()

    return args


class System(object):

    def __init__(self, initial, nsolutes=1, ff='gaff', T=300, pcoupltype='semiisotropic', mpi=False, np=4):
        """
        :param initial: initial configuration with no solutes
        :param nsolutes : number of solutes to add to each pore
        :param equil_length: length of equlibration simulation
        :param umbrella_length: length of umbrella sampling simulations
        :param ff: forcefield to use
        :param T: simulation temperature
        :param mpi: set to True, if you are running gromacs in parallel
        :param np: number of processes (if mpi = True)
        """

        print('Adding solutes to initial configuration...', end='', flush=True)
        p = subprocess.Popen(['place_solutes_pores.py', '-g', '%s' % initial, '-o', 'solutes.gro',
                              '-n', '%s' % nsolutes])  # creates solute.top
        p.wait()
        print('Done!')

        # create topology from new system
        self.top = SystemTopology('solutes.gro', ff=ff)
        self.top_name = 'topol.top'
        self.top.write_top(name=self.top_name)

        self.mdp = SimulationMdp(initial, T=T, p_coupling=pcoupltype)
        self.mdp.write_em_mdp(out='em')

        if mpi:
            self.gmx = 'mpirun -np %s gmx_mpi' % np
        else:
            self.gmx = 'gmx'

    def energy_minimize(self, name, out='em'):

        p = subprocess.Popen(['%s' % self.gmx, 'grompp', '-f', 'em.mdp', '-p', '%s' % self.top_name, '-c', '%s' % name,
                              '-o', '%s' % out])
        p.wait()

        run_string = "%s mdrun -v -deffnm %s" % (self.gmx, out)
        p2 = subprocess.Popen(run_string, shell=True)
        p2.wait()

    def npt_equilibration(self, conf, mdp_name='npt', out='npt', length=1000, barostat='berendsen'):
        """
        :param conf: initial configuration to be simulated
        :param mdp_name: name of .mdp file to be created
        :param out: name of .tpr file
        :param length: simulation length (ps)
        """

        self.mdp.length = length  # modify simulation length
        self.mdp.barostat = barostat  # modify barostat

        maxwarn = 0
        if barostat == 'Parrinello-Rahman':
            maxwarn = 1  # GROMACS may complain about generating velocities

        self.mdp.write_npt_mdp(out=mdp_name)
        p = subprocess.Popen(['gmx', 'grompp', '-f', '%s.mdp' % mdp_name, '-p', '%s' % self.top_name, '-c', '%s' % conf,
                              '-o', '%s' % out, '-maxwarn', '%d' % maxwarn])
        p.wait()

        run_string = "%s mdrun -v -deffnm %s" % (self.gmx, out)
        p2 = subprocess.Popen(run_string, shell=True)
        p2.wait()


if __name__ == "__main__":

    args = initialize()

    # initialize system by adding solutes to initial configuration and preparing topologies and .mdp file objects
    sys = System(args.gro)
    #sys.energy_minimize('solutes.gro', out='solutes_em')  # energy minimize solutes.gro
    sys.npt_equilibration('solutes_em', mdp_name='npt_initial_equil', out='npt_initial_equil', length=args.equil_length,
                          barostat='Parrinello-Rahman')  # run npt equilibration