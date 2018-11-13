#! /usr/bin/env python

import os
from builtins import range
from past.utils import old_div
import argparse
import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.analysis import Poly_fit, top, Atom_props, p2p
import mdtraj as md
import time
from scipy import stats


def initialize():

    parser = argparse.ArgumentParser(description='Calculate diffusion coefficient')
    parser.add_argument('-t', '--trajectory', default='wiggle.trr', help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of .gro coordinate file')
    parser.add_argument('-r', '--residue', type=str, help='Name of residue whose diffusivity we want')
    parser.add_argument('-atoms', nargs='+', help='Name of atoms whose collective diffusivity is desired')
    parser.add_argument('-b', '--nboot', default=200, help='Number of bootstrap trials to be run')
    parser.add_argument('-f', '--frontfrac', default=0.05, type=float, help='Where to start fitting line on msd curve')
    parser.add_argument('-F', '--fracshow', default=.2, type=float, help='Percent of graph to show, also where to stop '
                        'fitting line during diffusivity calculation')
    parser.add_argument('-a', '--axis', default='xyz', type=str, help='Which axis to compute msd along')
    parser.add_argument('-nboot', default=200, type=int, help='Number of bootstrap trials for error estimation')
    # | not implemented |
    # v                 v
    parser.add_argument('--restrict_to_pores', action="store_true", help='Only look at residue within pores of HII'
                                                                         'membrane')
    parser.add_argument('-radius', default=1, type=float, help='Radius of pores. Anything greater than this distance'
                        'from the pore center will not be included in calculation')
    args = parser.parse_args()

    return args


def autocorrFFT(x):

    N = len(x)
    F = np.fft.fft(x, n=2*N)  # 2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real   # now we have the autocorrelation in convention B
    n= N*np.ones(N)-np.arange(0, N)  # divide res(m) by (N-m)
    return old_div(res,n)  # this is the autocorrelation in convention A


def msd_fft(x, ndx):

    r = np.copy(x)
    r = r[:, ndx]
    N = len(r)
    D = np.square(r).sum(axis=1)
    D = np.append(D, 0)
    S2 = sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    Q = 2*D.sum()
    S1 = np.zeros(N)
    for m in range(N):
      Q = Q-D[m-1]-D[N-m]
      S1[m] = old_div(Q,(N-m))
    return S1-2*S2


def msd(x, axis):
    """
    Calculate mean square displacement based on particle positions
    :param x: particle positions
    :param axis: axis along which you want MSD (x, y, z, xy, xz, yz, xyz)
    :return: MSD of each particle
    """

    nT = x.shape[0]
    no_comp = x.shape[1]
    # MSD = np.zeros([nT], dtype=float)
    MSDs = np.zeros([nT, no_comp], dtype=float)  # a set of MSDs per particle

    for n in range(no_comp):
        MSDs[:, n] = msd_fft(x[:, n, :], axis)
    #     MSD += MSDs[:, n]
    # MSD /= no_comp

    return MSDs


def msd_straightforward(x, ndx):
    """
    Straightforward way to calculte msd. Gives same answer as msd()
    :param x: positions of centers of mass of all particles for each frame, numpy array [nframes, natoms, dim]
    :param ndx: list of indices to include in msd calculation (x = 0, y = 1, z = 2)
    :return: Average MSD and individual particle MSDs
    """
    n = x.shape[1]  # number of atoms
    N = x.shape[0]  # number of frames

    MSD = np.zeros([N])
    MSDs = np.zeros([N, n])

    for m in range(N):  # there nT different length intervals we can look at
        for k in range(N - m - 1):  # there are N - m - 1 independent intervals of length m
            MSDs[m, :] += np.linalg.norm(x[k + m, :, ndx] - x[k, :, ndx], axis=0)**2  # mean square displacement of all particles summed over all intervals of length m
        MSDs[m, :] /= (N - m)  # divide by the number of intervals to get an average msd over length m
        MSD[m] = np.mean(MSDs[m, :])

    return MSD, MSDs


class Diffusivity(object):

    def __init__(self, traj, gro, axis, begin=0, startfit=0.05, endfit=0.2, residue=False, atoms=[], restrict=[]):
        """
        Calculate diffusivity from trajectory
        :param traj: unwrapped trajectory (i.e. gmx trjconv with -pbc nojump)
        :param gro: representative coordinate file
        :param axis: axis along which to compute MSD
        :param startfit: start linear fit to MSD startfit % into trajectory
        :param endfit: end linear fit to MSD endfit % into trajectory
        :param residue: if specified, the residue whose center of mass MSD will be measured
        :param atoms: if specified, group of atoms whose center of mass MSD will be measured
        :param restrict: restrict selection to certain indices. For example, if you want to calculate MSD of a certain
        residue, but only include a fraction of the total residues in the system.
        """

        self.script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        self.top_location = "%s/../top/topologies" % self.script_location

        print('Loading trajectory...', end='', flush=True)
        t = md.load(traj, top=gro)[begin:]  # load trajectory
        print('Done!')
        self.nT = t.n_frames  # number of frames
        self.time = t.time  # time stamp on each frame
        self.MSD = None
        self.limits = None
        self.MSD_average = 0
        self.startfit = int(startfit*self.nT)  # index at which to start fit
        self.endfit = int(endfit*self.nT)  # index at which to end fit
        self.y_fit = []
        self.A = []  # fit parameters
        self.W = []  # weight matrix for curve fitting
        self.errorevery = int(np.ceil(self.nT / 100.0))  # plot only 100 bars total
        self.confidence_interval = 0
        self.yfit = None
        self.slope_error = 0
        self.Davg = 0

        if residue:

            res = residue
            if res == 'SOL':  # mdtraj changes the name from SOL to HOH
                res = 'HOH'

            self.itp = "%s/%s.itp" % (self.top_location, res)

            if restrict:
                selection = [a.index for a in t.topology.atoms if a.residue.name == res and a.index in restrict]
            else:
                selection = [a.index for a in t.topology.atoms if a.residue.name == res]

            topology = top.Top(self.itp)  # read topology
            atoms_per_residue = topology.natoms  # number atoms in a single residue
            matoms = topology.atom_masses  # mass of the atoms in residue
            self.mres = np.sum(matoms)  # total mass of residue

        elif atoms:

            selection = [a.index for a in t.topology.atoms if a.name in atoms]

            atoms_per_residue = len(atoms)
            matoms = np.array([Atom_props.mass[x] for x in atoms])
            self.mres = np.sum(matoms)
        else:
            print('Error: No valid group of atoms or residues selected')
            exit()

        pos = t.xyz[:, selection, :]

        self.axis = []
        if 'x' in axis:
            self.axis.append(0)
        if 'y' in axis:
            self.axis.append(1)
        if 'z' in axis:
            self.axis.append(2)

        print('Calculating center of mass of residues...', end='', flush=True)
        self.com = np.zeros([self.nT, pos.shape[1] // atoms_per_residue, 3])  # track the center of mass of each residue

        for f in range(self.nT):
            for i in range(self.com.shape[1]):
                w = (pos[f, i * atoms_per_residue:(i + 1) * atoms_per_residue, :].T * matoms).T  # weight each atom in the residue by its mass
                self.com[f, i, :] = np.sum(w, axis=0) / self.mres  # sum the coordinates and divide by the mass of the residue
        print('Done!')

        self.weights = True
        if self.com.shape[1] == 1:
            self.weights = False

        self.dt = t.time[-1] - t.time[-2]  # time step (assuming equispaced time points)

    def calculate(self):

        print('Calculating MSD...', end='', flush=True)
        self.MSD = msd(self.com, self.axis)
        self.MSD_average = np.mean(self.MSD, axis=1)
        print('Done!')

    def ensure_fit(self):

        fit = 0
        while fit == 0:

            self.yfit, _, self.slope_error, _, A = Poly_fit.poly_fit(self.time[self.startfit:self.endfit],
                                                  self.MSD_average[self.startfit:self.endfit], 1, self.W)

            plt.plot(self.time[self.startfit:self.endfit]/1000, self.yfit, '--', color='black', label='Linear Fit')

            # plt.errorbar(self.time, self.MSD_average, yerr=[self.limits[0, :], self.limits[1, :]],
            #              errorevery=self.errorevery, label='MSD')
            plt.plot(self.time/1000, self.MSD_average, label='MSD')

            plt.ylabel('MSD ($nm^2$)', fontsize=14)
            plt.xlabel('time (ns)', fontsize=14)
            plt.gcf().get_axes()[0].tick_params(labelsize=14)
            plt.legend(loc=2)
            plt.tight_layout()
            plt.ion()
            plt.show()
            fit = int(input("Type '1' if the fit looks good: "))
            if fit != 1:
                print('Press enter to following prompts to leave as is')
                self.startfit = float(input("Time to start fit (ns): ") or self.startfit)
                self.endfit = float(input("Time to stop fit (ns): ") or self.endfit)
                self.startfit = int(self.startfit / (self.dt/1000.0))  # convert time to index in t.time
                self.endfit = int(self.endfit / (self.dt/1000.0))
            plt.clf()

    def bootstrap(self, N):
        """
        Estimate error at each point in the MSD curve using bootstrapping
        :param N: number of bootstrap trials
        """

        eMSDs = np.zeros([self.nT, N], dtype=float)  # create n bootstrapped trajectories

        for b in range(N):
            indices = np.random.randint(0, self.com.shape[1], self.com.shape[1])  # randomly choose particles with replacement
            for n in range(self.com.shape[1]):
                eMSDs[:, b] += self.MSD[:, indices[n]]  # add the MSDs of a randomly selected particle
            eMSDs[:, b] /= self.com.shape[1]  # Divide every timestep by Nparticles -- average the MSDs

        self.limits = np.zeros([2, self.nT], dtype=float)  # upper and lower bounds at each point along MSD curve
        # determine a 95 % error bound for each tau (out of n MSD's, use that for the error bars)
        for t in range(self.nT):
            self.limits[0, t] = np.abs(np.percentile(eMSDs[t, :], 2.5) - self.MSD_average[t])
            self.limits[1, t] = np.abs(np.percentile(eMSDs[t, :], 97.5) - self.MSD_average[t])

        npts = self.endfit - self.startfit
        if self.weights:
            self.W = np.zeros((npts, npts))
            for i in range(npts):
                self.W[i, i] = old_div(1, ((self.limits[0, i + self.startfit]) ** 2))
        else:
            self.W = 'none'

        slopes = np.zeros([N])
        for b in range(N):
            # fit line to each bootstrapped MSD
            A = Poly_fit.poly_fit(self.time[self.startfit:self.endfit], eMSDs[self.startfit:self.endfit, b],
                                  1, self.W)[-1]
            slopes[b] = A[1]

        slopes /= (2*1000000*len(self.axis))

        # calculate 95 % confidence interval
        self.confidence_interval = stats.t.interval(0.95, N - 1, loc=slopes.mean(), scale=stats.sem(slopes))
        self.Davg = slopes.mean()

    def plot(self, axis, fracshow=0.5, savedata=False):

        end_frame = int(fracshow * self.time.size)
        #plt.plot(self.time[self.startfit:self.endfit]/1000, self.yfit, '--', color='black', label='Linear Fit')
        plt.errorbar(self.time[:end_frame]/1000, self.MSD_average[:end_frame], yerr=[self.limits[0, :end_frame],
                     self.limits[1, :end_frame]], errorevery=self.errorevery, label='MSD')

        if savedata:

            np.savez_compressed('msd.npz', time=self.time[:end_frame]/1000, msd=self.MSD_average[:end_frame],
                                yerr=[self.limits[0, :end_frame], self.limits[1, :end_frame]])

        plt.ylabel('MSD ($nm^2$)', fontsize=14)
        plt.xlabel('time (ns)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        #plt.title('D = %1.2e $\pm$ %1.2e $m^{2}/s$' % (self.Davg, np.abs(self.Davg - self.confidence_interval[0])))
        #plt.legend(loc=2)
        plt.tight_layout()
        plt.savefig('Diffusivity_%s.pdf' % axis)
        plt.show()


if __name__ == "__main__":

    args = initialize()

    D = Diffusivity(args.trajectory, args.gro, args.axis, residue=args.residue, atoms=args.atoms)
    D.calculate()
    D.ensure_fit()  # make sure diffusivity is being measured from the linear region of the MSD curve
    D.bootstrap(args.nboot)
    D.plot(args.axis, fracshow=args.fracshow)
    print('D = %1.2e +/- %1.2e m^2/s' % (D.Davg, np.abs(D.Davg - D.confidence_interval[0])))
