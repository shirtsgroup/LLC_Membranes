#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import sys
from LLC_Membranes.llclib import topology, physical, file_rw


def initialize():

    parser = argparse.ArgumentParser(description='Use geometric criteria to identify hydrogen bonds')

    parser.add_argument('-t', '--traj', default='PR_nojump.xtc', type=str, help='GROMACS trajectory file (.xtc or .trr)')
    parser.add_argument('-g', '--gro', default='berendsen.gro', type=str, help='GROMACS coordinate file (.gro)')
    parser.add_argument('-m', '--monomer', default='NAcarb11V', type=str, help='Name of liquid crystal monomer used to'
                                                                               'build unit cell (no extension)')
    parser.add_argument('-n', '--npts', default=10, type=int, help='Number of points in spline')
    parser.add_argument('-l', '--load', default=False, help='Load pickled system object')
    parser.add_argument('-s', '--savename', default='hbonds.pl', type=str, help='Name of pickled object to save after'
                                                                                'calcualtions are complete')

    return parser


class Spline(object):

    def __init__(self, gro, traj, monomer, progress=True, npts_spline=10):
        """ Calculate the positions of a spline running through each pore of an LLC Membrane

        :param gro: name of GROMACS coordinate file
        :param traj: name of GROMACS trajectory file
        :param monomer: name of liquid crystal monomer used to build system
        :param progress: show a progress bar while spline is built
        :param npts_spline: number of points making up the spline in each pore

        :type gro: str
        :type traj: str
        :type monomer: str
        :type progress: bool
        :type npts_spline: int
        """

        print('Loading trajectory...', flush=True, end='')
        self.t = md.load(traj, top=gro)
        print('Done!')

        self.monomer = topology.LC(monomer)

        pore_defining_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.pore_defining_atoms
                               and a.residue.name in self.monomer.residues]

        if not pore_defining_atoms:
            sys.exit('There are no atoms specified which can be used to define the pore centers. Did you specify '
                     'the correct monomer for this system? If so, have you properly annotated the pore defining '
                     'atoms?')

        self.pore_centers = physical.avg_pore_loc(4, self.t.xyz[:, pore_defining_atoms, :], self.t.unitcell_vectors,
                                                  spline=True, progress=progress, npts=npts_spline)

    def build_spline(self, rep='K', frame=-1, name='spline'):
        """ Build the spline into the last frame of the trajectory

        :param rep: name of atom to use to represent spline
        :param frame: Index of frame to draw spline through
        :param name: name of output file (will be in .gro format)

        :type rep: str
        :type frame: int
        :type name: str

        :return: 'spline.gro'
        """
        pos = self.t.xyz[frame, ...]

        for i in range(4):
            pos = np.concatenate((pos, self.pore_centers[frame, i, ...]))
        ids = [a.name for a in self.t.topology.atoms]
        res = [a.residue.name for a in self.t.topology.atoms]
        ids += [rep] * self.pore_centers.shape[2] * self.pore_centers.shape[1]
        res += [rep] * self.pore_centers.shape[2] * self.pore_centers.shape[1]

        file_rw.write_gro_pos(physical.wrap_box(pos, self.t.unitcell_vectors[frame, ...]), 'spline.gro',
                              ucell=self.t.unitcell_vectors[frame, ...], ids=ids, res=res)


if __name__ == "__main__":

    args = initialize().parse_args()

    spline = Spline(args.gro, args.traj, args.monomer, npts_spline=args.npts)
    spline.build_spline()
