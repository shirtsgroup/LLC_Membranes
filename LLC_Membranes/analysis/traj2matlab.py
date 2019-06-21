#!/usr/bin/env python

"""
Save solute trajectories into objects that can be loaded by MATLAB
"""

import argparse
import numpy as np
import mdtraj as md
import scipy.io as io
from LLC_Membranes.llclib import physical, topology


def initialize():

    parser = argparse.ArgumentParser(description='Calculate and plot slice of full 3D correlation function. Currently,'
                                                 ' only slices in the z direction for a monoclinic unit cell are '
                                                 ' implemented.')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Trajectory file. Make sure to '
                                                                            'preprocess with gmx trjconv -pbc whole')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-a', '--atoms', nargs='+', action='append', help='Name of atoms to calculate correlation '
                        'function with respect to. The center of mass will be used')
    parser.add_argument('-r', '--res', help='Residue to create correlation function with '
                        'respect to. Will use center of mass. Will override atoms')
    parser.add_argument('-b', '--begin', default=0, type=int, help='First frame to load')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Last frame to load')
    parser.add_argument('-s', '--skip', default=1, type=int, help='Skip every `skip` frames when loading trajectory')

    return parser


class Traj2Matlab(object):

    def __init__(self, traj, gro, begin=0, end=-1, skip=1):
        """ Initialize trajectory

        :param traj: gromacs trajectory (.xtc or .trr)
        :param gro: gromacs coordinate file (.gro)
        :param begin: first frame to load
        :param end: last frame to load
        :param skip: skip every `skip` frames

        :type traj: str
        :type gro: str
        :type begin: int
        :type end: int
        :type skip: int
        """

        print('Loading Trajectory...', end='', flush=True)
        if end == -1:
            self.t = md.load(traj, top=gro)[begin::skip]
        else:
            self.t = md.load(traj, top=gro)[begin:end:skip]  # this excludes last frame
        print('Done!')

        self.residue_trajectory = None

    def extract_residue_trajectory(self, resname):
        """ Extract the center of mass coordinates of a specific residue

        :param resname: name of residue whose trajectory is desired

        :type resname: str
        """

        residue = topology.Residue(resname)
        resatoms = [a.index for a in self.t.topology.atoms if a.residue.name == resname]
        atom_names = [a.name for a in self.t.topology.atoms if a.residue.name == resname][:residue.natoms]
        mass = [residue.mass[a] for a in atom_names]
        self.residue_trajectory = physical.center_of_mass(self.t.xyz[:, resatoms, :], mass)

    def write(self, name='trajectory.mat'):
        """ Write trajectory to MATLAB-readable file (.mat)

        :param name: name of file to write

        :type name: str
        """

        io.savemat(name, dict(traj=self.residue_trajectory))


if __name__ == "__main__":

    args = initialize().parse_args()

    traj = Traj2Matlab(args.traj, args.gro, begin=args.begin, end=args.end, skip=args.skip)
    traj.extract_residue_trajectory(args.res)
    traj.write('trajectory_%s.mat' % args.res)
