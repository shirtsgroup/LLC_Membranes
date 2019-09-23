#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import sys
from LLC_Membranes.llclib import topology, physical, file_rw
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Use geometric criteria to identify hydrogen bonds')

    parser.add_argument('-t', '--traj', default='PR_nojump.xtc', type=str, help='GROMACS trajectory file (.xtc or .trr)')
    parser.add_argument('-g', '--gro', default='berendsen.gro', type=str, help='GROMACS coordinate file (.gro)')
    parser.add_argument('-m', '--monomer', default='NAcarb11V', type=str, help='Name of liquid crystal monomer used to'
                                                                               'build unit cell (no extension)')
    parser.add_argument('-n', '--npts', default=10, type=int, help='Number of points in spline')
    parser.add_argument('-l', '--load', default=False, help='Load pickled system object')
    parser.add_argument('-s', '--savename', default='spline.pl', type=str, help='Name of pickled object to save after'
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

        self.tortuosity = None
        self.length = None

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

    def compute_tortuosity(self):
        """ Calculate the tortuosity of the pores by computing the ratio of L / Z where L is the length of the pore
        spline curve and Z is the length of the unit cell in the z direction

        """

        npts = self.pore_centers.shape[2]
        npores = self.pore_centers.shape[1]
        self.tortuosity = np.zeros([self.t.n_frames, npores])
        self.length = np.zeros([self.t.n_frames, npores])

        for t in range(self.t.n_frames):

            zbox = self.t.unitcell_vectors[t, 2, 2]  # z box dimension for this frame

            l = np.zeros(npores)
            for p in range(npts - 1):
                l += np.linalg.norm(self.pore_centers[t, :, p + 1, :] - self.pore_centers[t, :, p, :], axis=1)

            # Minimum image distance between 1st and last spline points
            bot = self.pore_centers[t, :, 0, :] + [0, 0, zbox]  # shift bottom point up
            l += np.linalg.norm(bot - self.pore_centers[t, :, -1, :], axis=1)
            self.tortuosity[t, :] = l / zbox
            self.length[t, :] = l

    def plot_tortuosity(self):

        npores = self.pore_centers.shape[1]

        for i in range(npores):
            plt.plot(self.t.time, self.tortuosity[:, i], linewidth=2)

        plt.show()

    def save_tortuosity(self):

        file_rw.save_object(self.tortuosity, 'tortuosity.pl')


if __name__ == "__main__":

    args = initialize().parse_args()

    spline = Spline(args.gro, args.traj, args.monomer, npts_spline=args.npts)
    spline.build_spline()
    spline.compute_tortuosity()
    print('Mean Tortuosity: %.2f +/- %.2f' % (spline.tortuosity.mean(), spline.tortuosity.std()))
    print('Mean Pore Length: %.2f +/- %.2f' % (spline.length.mean(), spline.length.std()))
    spline.save_tortuosity()

