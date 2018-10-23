#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import Atom_props
import ruptures
from scipy.optimize import curve_fit
from LLC_Membranes.analysis.hbonds import Residue


def initialize():

    parser = argparse.ArgumentParser(description='Model a continuous time random walk')

    # MD trajectory control
    parser.add_argument('-t', '--trajectory', default='PR_nojump.xtc', help='Path to input file.')
    parser.add_argument('-g', '--gro', default='em.gro', help='Name of .gro coordinate file.')
    parser.add_argument('-b', '--begin', default=0, help='First trajectory frame used for analysis.')
    parser.add_argument('-e', '--end', default=-1, help='Last trajectory frame used for analysis.')
    parser.add_argument('-step', '--step', default=1, help='Include every "step" frames')

    parser.add_argument('-r', '--residue', default='ETH', type=str, help='Name of residue whose diffusivity we want')
    parser.add_argument('-atoms', nargs='+', help='Name of atoms whose collective diffusivity is desired')
    parser.add_argument('-nbins', default=25, type=int, help='Number of bins to bin hop and dwell distributions into')

    # ctrw simulation
    parser.add_argument('-n', '--nhops', default=10000, type=int, help='Number of hops to perform')

    # bootstrapping
    # parser.add_argument('-b', '--nboot', default=200, help='Number of bootstrap trials to be run')
    # parser.add_argument('-f', '--frontfrac', default=0.05, type=float, help='Where to start fitting line on msd curve')
    # parser.add_argument('-F', '--fracshow', default=.2, type=float, help='Percent of graph to show, also where to stop '
    #                     'fitting line during diffusivity calculation')
    # parser.add_argument('-a', '--axis', default='xyz', type=str, help='Which axis to compute msd along')
    # parser.add_argument('-nboot', default=200, type=int, help='Number of bootstrap trials for error estimation')
    # # | not implemented |
    # # v                 v
    # parser.add_argument('--restrict_to_pores', action="store_true", help='Only look at residue within pores of HII'
    #                                                                      'membrane')
    # parser.add_argument('-radius', default=1, type=float, help='Radius of pores. Anything greater than this distance'
    #                     'from the pore center will not be included in calculation')

    return parser


def autocorrFFT(x):

    N = len(x)
    F = np.fft.fft(x, n=2*N)  # 2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real   # now we have the autocorrelation in convention B
    n= N*np.ones(N)-np.arange(0, N)  # divide res(m) by (N-m)

    return res / n  # this is the autocorrelation in convention A


def msd_fft(x):

    N = len(x)
    D = np.square(x)
    D = np.append(D, 0)
    S2 = autocorrFFT(x)
    Q = 2*D.sum()
    S1 = np.zeros(N)
    for m in range(N):
      Q = Q-D[m-1]-D[N-m]
      S1[m] = Q / (N-m)

    return S1-2*S2


def gaussian(points, mean, sigma, amplitude):
    return (amplitude / np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(
        -(points - mean) ** 2 / (2 * sigma ** 2))


def wait_times(t, A, lam):

    return A*np.exp(-lam * t)


def wait_times_empirical(t, A):

    return A*np.exp(-lambda_empirical * t)


def wait_times_cdf(t, lam):

    return 1 - np.exp(-lam * t)


def random_dwell_time(lam):

        return -np.log(1 - np.random.uniform()) / lam


class System(object):

    def __init__(self, traj, gro, res, start=0, end=-1, step=1):

        self.t = md.load(traj, top=gro)[start:end:step]
        self.t.time /= 1000  # convert to nanoseconds
        self.dt = self.t.time[1] - self.t.time[0]

        keep = [a.index for a in self.t.topology.atoms if a.residue.name == res]

        self.residue = Residue(res)
        self.nres = len(keep) // self.residue.natoms
        self.mass = np.array([Atom_props.mass[a.name] for a in self.t.topology.atoms if a.residue.name == res][:self.residue.natoms])
        self.pos = self.t.xyz[:, keep, :]

        self.com = np.zeros([self.t.n_frames, self.nres, 3])

        for f in range(self.t.n_frames):
            for i in range(self.nres):
                w = (self.pos[f, (i * self.residue.natoms):((i + 1)*self.residue.natoms), :].T * self.mass).T
                self.com[f, i, :] = np.sum(w, axis=0) / sum(self.mass)

        # initialize for later
        self.dwell_times = []
        self.hop_lengths = []
        self.lambda_distribution = []  # distribution of lambda for poisson process
        self.hop_sigma_distribution = []  # distribution of standard deviation of hop lengths

    def hops_and_dwells(self, penalty=1):
        """ Find breakpoints then assemble lists of dwell times and hop lengths. See documentation for Ruptures python
        package: http://ctruong.perso.math.cnrs.fr/ruptures-docs/build/html/index.html for more options that can be
        added

        :param penalty: penalty for cost function minimization
        :return:
        """

        for j in range(self.nres):

            bp = ruptures.detection.Pelt().fit_predict(self.com[:, j, :], pen=penalty)

            for i in range(len(bp) - 2):  # exclude first and last segments
                self.dwell_times.append(self.t.time[bp[i + 1]] - self.t.time[bp[i]])
                if i > 0:
                    self.hop_lengths.append(np.mean(self.com[bp[i]:bp[i + 1], j, 2]) -
                                            np.mean(self.com[bp[i - 1]:bp[i], j, 2]))

    def fit_distributions(self, nbins=25, plot=True, nboot=200):

        # Reconstruct bootstrapped distributions by randomly sampling from data
        # keep all bootstrapping data to generate error bars
        hist_dwell = np.zeros([nboot, nbins])
        bins_dwell_centered = np.zeros([nboot, nbins])
        p_dwell = [1, len(self.dwell_times) / (self.nres * self.t.time[-1])]
        hist_jump = np.zeros([nboot, nbins])
        bins_jump_centered = np.zeros([nboot, nbins])

        for i in range(nboot):

            # dwell times
            dwell_times_boot = np.random.choice(self.dwell_times, size=len(self.dwell_times), replace=True)
            hist_dwell[i, :], bins_dwell = np.histogram(dwell_times_boot, bins=nbins)
            bins_dwell_centered[i, :] = np.array([(bins_dwell[i + 1] + bins_dwell[i])/2 for i in
                                                  range(len(hist_dwell[i, :]))])
            solp_dwell, cov_x = curve_fit(wait_times, bins_dwell_centered[i, :], hist_dwell[i, :], p_dwell)

            self.lambda_distribution.append(solp_dwell[1])

            # jump lengths
            hop_lengths_boot = np.random.choice(self.hop_lengths, size=len(self.hop_lengths), replace=True)
            hist_jump[i, :], bins_jump = np.histogram(hop_lengths_boot, bins=nbins)
            bins_jump_centered[i, :] = np.array([(bins_jump[i + 1] + bins_jump[i])/2 for i in
                                                 range(len(hist_jump[i, :]))])
            p_hops = [np.mean(hist_jump[i, :]), np.std(hist_jump[i, :]), 1]
            solp_hops, cov_x = curve_fit(gaussian, bins_jump_centered[i, :], hist_jump[i, :], p_hops)

            self.hop_sigma_distribution.append(np.abs(solp_hops[1]))  # sometime curve_fit finds negative answer

        if plot:
            # change so that solp's used below are from average parameters
            fig, ax = plt.subplots(2, 2, figsize=(12, 8))

            ax[0, 0].bar(bins_dwell_centered.mean(axis=0), hist_dwell.mean(axis=0),
                         width=(bins_dwell[1] - bins_dwell[0]))
            ax[0, 0].set_xlabel('Dwell Time (ns)')
            ax[0, 0].set_ylabel('Frequency')
            ax[0, 0].plot(bins_dwell_centered.mean(axis=0), wait_times(bins_dwell_centered.mean(axis=0), solp_dwell[0],
                       solp_dwell[1]), '--', color='black', label='$\lambda_{fit}$ = %.3f' % solp_dwell[1])
            ax[0, 0].legend()

            ax[0, 1].hist(self.lambda_distribution, bins=nbins)

            ax[1, 0].bar(bins_jump_centered.mean(axis=0), hist_jump.mean(axis=0), width=(bins_jump[1] - bins_jump[0]))
            ax[1, 0].plot(bins_jump_centered.mean(axis=0), gaussian(bins_jump_centered.mean(axis=0), solp_hops[0],
                       solp_hops[1], solp_hops[2]), '--', label='Gaussian fit\n $\mu$=%.2f, $\sigma$=%.2f' %
                                                                (solp_hops[0], solp_hops[1]), color='black')
            ax[1, 0].set_xlabel('Hop Length ($z$-direction, nm)')
            ax[1, 0].set_ylabel('Frequency')
            ax[1, 0].legend()

            ax[1, 1].hist(self.hop_sigma_distribution, bins=nbins)


class CTRW(object):

    def __init__(self, nhops, dwell_lambdas, hop_sigmas, ntrials):

        self.dwell_lambdas = dwell_lambdas
        self.hop_sigmas = hop_sigmas
        self.nhops = nhops

        self.trajectories = np.zeros([ntrials, nhops, 2])  # last dimension is [time, z_position]
        self.ntrials = ntrials
        self.trial = 0

        self.trajectory_hops = np.zeros([ntrials, 2*self.nhops - 1, 2])


    def generate_trajectory(self):

        # constrain mean to be zero
        z_position = np.cumsum(np.random.normal(loc=0, scale=np.random.choice(self.hop_sigmas), size=self.nhops))
        self.trajectories[self.trial, :, 1] = z_position - z_position[0]  # make initial z equal to 0

        time = np.zeros([self.nhops])
        lambda_trial = np.random.choice(self.dwell_lambdas)
        for i in range(1, self.nhops):  # make initial time equal to 0
            time[i] = random_dwell_time(lambda_trial)  # hop at random time intervals according to poisson process
        self.trajectories[self.trial, :, 0] = np.cumsum(time)

        self.trajectory_hops[self.trial, 1::2, 0] = time[1:]
        self.trajectory_hops[self.trial, 2::2, 0] = time[1:]

        self.trajectory_hops[self.trial, ::2, 1] = z_position
        self.trajectory_hops[self.trial, 1:-1:2, 1] = z_position[:-1]
        self.trajectory_hops[self.trial, -1, 1] = z_position[-1]


if __name__ == "__main__":

    args = initialize().parse_args()

    sys = System(args.trajectory, args.gro, args.residue, start=args.begin, end=args.end, step=args.step)
    sys.hops_and_dwells()
    sys.fit_distributions(nbins=args.nbins)

    random_walks = CTRW(args.nhops, sys.lambda_distribution, sys.hop_sigma_distribution, 200)

    # make uniform time intervals
    n_pts = int(0.1*nhops)
    time_uniform = np.linspace(0, time[-1], n_pts)
    z_interpolated = np.zeros([n_pts, 1])
    for i, x in enumerate(time_uniform):
        time_index = np.argmin(np.abs(x - time))
        if x - time[time_index] < 0:
            time_index -= 1
        z_interpolated[i, 0] = z_position[time_index]
    #	print(x, time[time_index], z_interpolated[i])


    msd = msd_fft(z_interpolated[:, 0])

    fig, ax = plt.subplots(1, 2, figsize=(12,5))
    ax[0].plot(time, z_position, linewidth=2, label='Raw generated data')
    ax[0].plot(time_plottable, z_positions_plottable, '--', color='black', linewidth=2, label='Data converted to hops')
    ax[0].plot(time_uniform, z_interpolated[:, 0], label='Uniform time intervals')
    ax[0].set_xlabel('Time (ns)', fontsize=16)
    ax[0].set_ylabel('$z$-coordinate (nm)', fontsize=16)
    ax[0].legend()

    ax[1].plot(time_uniform, msd)
    #ax[1].plot(time_uniform, msd2)
    plt.show()
