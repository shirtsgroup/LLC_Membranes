#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import ruptures
from scipy.optimize import curve_fit, minimize
from LLC_Membranes.llclib import topology, physical, file_rw, fitting_functions, timeseries
import tqdm
from ctrwsim import CTRW
import sqlite3 as sql
import os


def initialize():

    parser = argparse.ArgumentParser(description='Model a continuous time random walk')

    # MD trajectory control
    parser.add_argument('-t', '--trajectory', default='PR_nojump.xtc', help='Path to input file.')
    parser.add_argument('-g', '--gro', default='em.gro', help='Name of .gro coordinate file.')
    parser.add_argument('-b', '--begin', default=0, help='First trajectory frame used for analysis.')
    parser.add_argument('-e', '--end', default=-1, help='Last trajectory frame used for analysis.')
    parser.add_argument('-step', '--step', default=1, help='Include every "step" frames')
    parser.add_argument('-ma', '--moving_average', default=False, type=int, help='Calculate a moving average of the '
                        'center of mass coordinate')

    # loading saved objects
    parser.add_argument('-load', '--load', default=False, help='Specify name of pickled object to load')

    # restrict to pores parameters
    parser.add_argument('-restrict', '--restrict_to_pores', action='store_true', help='Only look at residues which'
                        'stay in the pore (based on last simulation frame)')
    parser.add_argument('--pore_defining_atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Atoms'
                        'used to define pore centers')
    parser.add_argument('--pore_defining_residue', default='HII', type=str, help='residue to which pore_defining atoms '
                                                                                 'belong')
    parser.add_argument('-radius', '--pore_radius', default=1.48, type=float, help='Defined pore radius (nm)')

    parser.add_argument('-r', '--residue', default='ETH', type=str, help='Name of residue whose diffusivity we want')
    parser.add_argument('-atoms', nargs='+', help='Name of atoms whose collective diffusivity is desired')
    parser.add_argument('-nbins', default=25, type=int, help='Number of bins to bin hop and dwell distributions into')
    parser.add_argument('-bp', '--breakpoint_penalty', default=0.25, type=float, help='Cost function penalty when '
                                                                                      'determing break points')

    # ctrw simulation
    parser.add_argument('-ntsim', '--ntrajsim', default=1000, type=int, help='Number of trajectories to simulate')
    parser.add_argument('--update', action="store_true", help="update database with all parameters calculated")
    parser.add_argument('--ensemble', action="store_true", help="Calculate ensemble average MSD. If false, will do"
                                                                "time-averaged MSD")

    # bootstrapping
    parser.add_argument('-nboot', '--nboot', default=200, type=int, help='Number of bootstrap trials to be run')
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


def gaussian(points, mean, sigma):

    return (1 / np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(
        -(points - mean) ** 2 / (2 * sigma ** 2))


class System(object):

    def __init__(self, traj, gro, res, start=0, end=-1, step=1, ma=False):
        """ Using an MD trajectory, calculate fit parameters for the dwell time distribution and hop length distribution

        :param traj: name of MD trajectory (GROMACS .xtc or .trr)
        :param gro: name of .gro (or .pdb) coordinate file associated with same topology as traj
        :param res: name of residue that is hopping around
        :param start: first frame of trajectory to include
        :param end: last frame of trajectory to include
        :param step: include every step frames in calculations
        :param ma: calculate a moving average on the center of mass coordinate traces

        :type traj: str
        :type gro: str
        :type res: str
        :type start: int
        :type end: int
        :type step: int
        """

        # load trjaectory
        print('Loading trajectory...', end='', flush=True)
        self.t = md.load(traj, top=gro)[start:end:step]
        print('Done!')

        self.time = self.t.time / 1000  # time in nanoseconds
        self.dt = self.time[1] - self.time[0]  # time step

        keep = [a.index for a in self.t.topology.atoms if a.residue.name == res]  # get indices of atoms making up residue of interest

        self.res_start = keep[0]  # index where residue starts

        self.residue = topology.Residue(res)
        self.nres = len(keep) // self.residue.natoms  # number of residues in system
        self.mass = [v for v in self.residue.mass.values()]  # mass of atoms making up residue
        self.pos = self.t.xyz[:, keep, :]  # positions of residue atoms

        print('Calculating centers of mass...', end='', flush=True)
        self.com = physical.center_of_mass(self.pos, self.mass)  # center of mass of residues
        print('Done!')

        if ma:
            self.calculate_moving_average(ma)

        # initialize for later
        self.pore_centers = None
        self.dwell_times = []  # distribution of dwell times
        self.tail_dwells = []
        self.hop_lengths = []
        self.hurst_distribution = []
        self.hop_acf = None
        self.alpha_distribution = []  # distribution of alpha for poisson process
        self.hop_sigma_distribution = []  # distribution of standard deviation of hop lengths
        self.breakpoint_penalty = 0
        self.location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # This script location
        self.r = []  # radial distance of solute from pore center when hops are made

        self.partition = None  # array telling whether a solute is in the pores or tails. True if in pores, else False

    def plot_random_trajectories(self, n=1):
        """ Plot random center of mass trajectories and print indices of atoms in residue defining that COM

        :param n: Number of random trajectories to plot

        :type n: int
        """

        for i in np.random.randint(self.nres, size=n):
            start = self.res_start + len(self.mass)*i
            end = start + len(self.mass)
            print('Inidices %s through %s' % (start, end))
            plt.plot(self.com[:, i, 2], linewidth=2)
            plt.xlabel('Frame')
            plt.ylabel('Coordinate')
            plt.show()

    def calculate_moving_average(self, n):
        """ Calculate moving average of a time series

        :param n: Number of previous points to average

        :type n: int
        """

        nT = self.com.shape[0]
        ma = np.zeros([nT - n + 1, self.com.shape[1], 3])

        for i in range(self.com.shape[1]):
            for d in range(3):
                ret = np.cumsum(self.com[:, i, d], dtype=float)
                ret[n:] = ret[n:] - ret[:-n]
                ma[:, i, d] = ret[n - 1:] / n

        self.com = ma

    def calculate_solute_partition(self, r=1.5, spline=False, membrane_residue='HII', write_tcl=True):
        """ Determine whether each COM is in the tail or pore region as a function of time

        :param r: radial distance from pore center where pore region transitions to tail region
        :param spline: calculate pore centers as function of z using a spline in each pore
        :param membrane_residue: if using spline, give the name of the liquid crystal residue used to make the membrane

        :type r: float
        :type spline: bool
        :type membrane_residue: str
        """

        membrane = topology.LC('%s.gro' % membrane_residue)  # object w/ attributes of LC making up membrane
        keep = [a.index for a in self.t.topology.atoms if a.name in membrane.pore_defining_atoms and a.residue.name
                == membrane.name]
        self.pore_centers = physical.avg_pore_loc(4, self.t.xyz[:, keep, :], self.t.unitcell_vectors, buffer=0,
                                                  spline=spline, npts=10, progress=True, bins=False)

        self.partition = physical.partition(self.com, self.pore_centers, r, unitcell=self.t.unitcell_vectors,
                                            spline=spline)

        if write_tcl:
            self.write_tcl()

    def write_tcl(self, frame=-1, name='partition.tcl'):
        """ Write a .tcl script to view solute partition in VMD
        """

        pore_indices = [self.res_start + self.residue.natoms * i for i in np.where(self.partition[frame, :])[0]]
        tail_indices = [self.res_start + self.residue.natoms * i for i in np.where(self.partition[frame, :] == False)[0]]

        with open(name, 'w') as f:
            f.write('color Display Background white\n')
            f.write('mol addrep 0\n')
            f.write('mol modselect 0 0 index')
            for i in pore_indices:
                f.write(' %s' % i)
            f.write('\n')
            f.write('mol modcolor 0 0 ColorID 0\n')
            f.write('mol modstyle 0 0 CPK 2.0 0.3 12.0 12.0\n')
            f.write('mol addrep 0\n')
            f.write('mol modselect 1 0 index')
            for i in tail_indices:
                f.write(' %s' % i)
            f.write('\n')
            f.write('mol modstyle 1 0 CPK 2.0 0.3 12.0 12.0\n')
            f.write('mol modcolor 1 0 ColorID 1\n')

    def write_spline_coordinates(self, frame=-1):
        """ write spline coordinates to a .gro file
        """

        coordinates = np.reshape(self.pore_centers[frame, ...], (self.pore_centers.shape[1]*self.pore_centers.shape[2],
                                                                 3))
        file_rw.write_gro_pos(coordinates, 'spline.gro', name='K')

    def hops_and_dwells(self, penalty=0.25, nframes_dwell=10, locations=False):
        """ Find breakpoints then assemble lists of dwell times and hop lengths. See documentation for Ruptures python
        package: http://ctruong.perso.math.cnrs.fr/ruptures-docs/build/html/index.html for more options that can be
        added

        :param penalty: penalty for cost function minimization
        :param nframes_dwell: number of frames that a solute needs to stay in the region of interest in order to be
        analyzed for hops.
        :param locations: if True, this function will return the distance from the pore center when hops are made

        :type penalty: float
        :type nframes_dwell: int
        :type location: bool

        """

        self.breakpoint_penalty = penalty
        print('penalty = %s' % self.breakpoint_penalty)

        # Handle initial and final dwell times since they don't have a true beginning or end. If the dwell times are
        # longer than some limit, then they will be included so that we don't miss out on dwells that last on the order
        # of the entire length of the trajectory
        initial_dwells = []  # length of time particle initially dwells before its first hop
        final_dwells = []  # length of time particle dwells up until end of trajectory

        for j in tqdm.tqdm(range(self.com.shape[1])):

            hop_lengths = []  # hop lengths for this solute

            # Get indices where res jumps between regions
            # See https://stackoverflow.com/questions/36894822/how-do-i-identify-sequences-of-values-in-a-boolean-array
            switch_points = np.argwhere(np.diff(self.partition[:, j])).squeeze().tolist()
            try:
                switch_points.append(self.partition.shape[0])
            except AttributeError:
                switch_points = [switch_points]
                switch_points.append(self.partition.shape[0])

            # Analyze sub-trajectories where solute is in defined pore region
            begin = 0
            for i, end in enumerate(switch_points[::2]):  # only at switch points after which solute is in pores
                if (end - begin) > nframes_dwell:  # only include hops/dwells from trajectories long enough to analyze
                    # np.savez_compressed('MET_trace.npz', z=self.com[begin:end, j, :])
                    # exit()
                    # movement_3d = np.linalg.norm(self.com[begin:end, j, :], axis=1)  # magnitude of distance travelled
                    bp = ruptures.detection.Binseg(model='l2', min_size=1, jump=1).fit_predict(self.com[begin:end, j, :]
                                                   , pen=self.breakpoint_penalty)
                    # ruptures.display(self.com[begin:end, j, :], bp, figsize=(10, 6))
                    # plt.show()
                    # exit()

                    if locations:

                        pore = j % self.pore_centers.shape[1]

                        if len(self.pore_centers.shape) == 4:  # spline

                            pore_centers = self.pore_centers[begin:end, pore, ...][bp[:-1]]
                            coms = self.com[begin:end, j, :][bp[:-1]]
                            box = self.t.unitcell_vectors[begin:end, ...][bp[:-1]]

                            for k in range(len(bp[:-1])):
                                self.r += physical.radial_distance_spline(pore_centers[k, ...], coms[k, :][np.newaxis, :],
                                                                          box[k]).tolist()

                        else:  # no spline

                            pore_centers = self.pore_centers[begin:end, pore, :][bp[:-1]]
                            coms = self.com[begin:end, j, :2][bp[:-1]]
                            self.r += np.linalg.norm(pore_centers - coms, axis=1).tolist()

                    # for visualization predicted z-hops on top of time series
                    # traj_hops = np.zeros([2*len(bp) - 1, 2])
                    #
                    # traj_hops[1::2, 0] = self.time[bp[:-1]]*2
                    # traj_hops[2::2, 0] = self.time[bp[:-1]]*2

                    initial_dwells.append(bp[0])
                    if len(bp) > 1:  # handle case where the only break point is at the end of the simulation
                        final_dwells.append(bp[-1] - bp[-2])

                    for k in range(len(bp) - 1):

                        if (k + 1) < (len(bp) - 1):  # exclude first and last segments for dwell times
                            self.dwell_times.append(bp[k + 1] - bp[k])  # discrete form

                        if k > 0:  # need two different z coordinates to get a hop length
                            hop_lengths.append(np.mean(self.com[bp[k]:bp[k + 1], j, 2]) -
                                                    np.mean(self.com[bp[k - 1]:bp[k], j, 2]))
                        else:
                            # bp[0] != 0, rather it is the first break point location. bp[-1] is the last point in the
                            # time series though
                            hop_lengths.append(np.mean(self.com[bp[k]:bp[k + 1], j, 2]) -
                                               np.mean(self.com[0:bp[k], j, 2]))  # hop at first segment

                        # traj_hops[2*k, 0] = bp[k]
                        # traj_hops[2*k + 1, 0] = bp[k + 1]
                        # traj_hops[2*k:(2*k + 2), 1] = np.mean(self.com[bp[k]:bp[k + 1], j, 2])

                    # ruptures.display(self.com[begin:end, j, :], bp, figsize=(10, 6))
                    # plt.plot(traj_hops[:, 0], traj_hops[:, 1], '--', color='black')
                    # plt.show()
                    # exit()
                try:
                    begin = switch_points[2*i + 1]
                except IndexError:
                    pass

            self.hop_lengths.append(hop_lengths)

        self.tail_dwells = np.array(initial_dwells + final_dwells)

    def fit_distributions(self, nbins=25, plot=True, nboot=200, show=False, save=True):
        """ Fit curves to dwell time and hop length distributions

        :param nbins: number of bins in histograms
        :param plot: make plots of the averaged distributions and bootstrapped parameter estimates
        :param nboot: number of bootstrap trials
        :param show: show plot of bootstrapped distributions
        :param save: save plot of bootstrapped distributions

        :type nbins: int
        :type plot: bool
        :type nboot: int
        :type show: bool
        :type save: bool
        """

        # Reconstruct bootstrapped distributions by randomly sampling from data
        # keep all bootstrapping data to generate error bars
        hist_dwell = np.zeros([nboot, nbins])
        hist_jump = np.zeros([nboot, nbins])

        all_hops = []
        for i in self.hop_lengths:
            all_hops += i

        hop_mean = 0

        # Define bins so histograms can be added together
        bins_dwell = np.linspace(np.min(self.dwell_times), np.amax(list(self.dwell_times) + list(self.tail_dwells)),
                                 nbins + 1)
        dwell_bin_width = bins_dwell[1] - bins_dwell[0]
        bins_dwell_centered = np.array([i + dwell_bin_width/2 for i in bins_dwell[:-1]])
        bins_hop = np.linspace(np.amin(all_hops), np.amax(all_hops), nbins + 1)
        hop_bin_width = bins_hop[1] - bins_hop[0]
        bins_hop_centered = np.array([i + hop_bin_width/2 for i in bins_hop[:-1]])

        # plot/save initial distribution
        # np.savez_compressed('dwells.npz', dwell_times=self.dwell_times)
        # hist, edges = np.histogram(self.dwell_times, range=(1, 25), bins=25)
        # normalized = hist / len(self.dwell_times)
        # plt.bar(edges[:-1], normalized, 1, align='edge')
        #
        # plt.figure()
        # plt.hist(self.dwell_times, bins=50)
        # plt.show()
        # exit()

        i = 0
        while i < nboot:

            # dwell times
            try:
                print('Bootstrapping ... Trial %s\r' % (i + 1), end='', flush=True)

                # Recreate dwell time distribution by randomly choosing from all dwell times with replacement
                dwell_times_boot = np.random.choice(self.dwell_times, size=len(self.dwell_times), replace=True)

                # Add long dwell times randomly selected from beginning and end of trajectory
                tail_dwells = self.tail_dwells[np.where(self.tail_dwells > max(dwell_times_boot))[0]]

                # add new tail dwells if there are any
                try:
                    dwell_times_boot = np.concatenate((np.random.choice(tail_dwells, size=len(tail_dwells),
                                                                        replace=True), dwell_times_boot))
                except ValueError:
                    pass

                # Maximum likelihood estimate of alpha
                args = (dwell_times_boot, int(min(dwell_times_boot)), True)  # (dwell times, xmin, maximize = True)
                maximum = minimize(fitting_functions.power_law_discrete_log_likelihood, 1.5, args=args,
                                   bounds=[(1.01, 3)]).x[0]

                self.alpha_distribution.append(maximum - 1)

                # Graphical representation of log-likelihood maximization
                # ll = []  # list of log-likelihoods for each alpha tested
                # alphas = np.linspace(1.01, 3, 100)
                # for a in alphas:
                #     ll.append(fitting_functions.power_law_discrete_log_likelihood(a, dwell_times_boot, 1))
                #
                # plt.plot(alphas, ll)
                # plt.show()

                if plot:

                    # bin the bins
                    hist_dwell[i, :], _ = np.histogram(dwell_times_boot, bins_dwell, normed=True)

                # Fit a line to log-log plot of dwell time distribution and compare it to MLE plot
                # _, alpha = fitting_functions.fit_power_law(np.array(bins_dwell_centered), hist_dwell[i, :])
                #
                # plt.plot(bins_dwell_centered, bins_dwell_centered ** alpha, color='blue', linewidth=3,
                #          label=r'$\alpha$=%.2f (LS)' % -alpha)
                #
                # # compare MLE to Least squares estimate
                # plt.bar(bins_dwell_centered, hist_dwell[i, :], dwell_bin_width, align='center')
                # plt.plot(bins_dwell_centered, bins_dwell_centered ** - maximum, color='red',
                #          linewidth=3, label=r'$\alpha$=%.2f (MLE)' % maximum)
                # plt.legend(fontsize=14)
                # plt.xlabel('Dwell Time (Frames)', fontsize=14)
                # plt.ylabel('Probability', fontsize=14)
                # plt.gcf().get_axes()[0].tick_params(labelsize=14)
                # plt.tight_layout()
                # plt.show()
                # exit()

                # jump lengths
                hop_lengths_boot = np.random.choice(all_hops, size=len(all_hops), replace=True)

                if plot:

                    hist_jump[i, :], _ = np.histogram(hop_lengths_boot, bins_hop, normed=True)

                # Maximum likelihood estimate of mean and sigma of hop length distribution
                # args = (hop_lengths_boot, True)  # (hop_lengths times, maximize = True)
                # guess = np.array([np.mean(hop_lengths_boot), np.std(hop_lengths_boot)])
                # hop_mle = minimize(fitting_functions.gaussian_log_likelihood, guess, args=args,
                #                    bounds=[[-np.inf, np.inf], [0, np.inf]]).x

                hop_mle = [np.mean(hop_lengths_boot), np.std(hop_lengths_boot)]  # mean and std are mle estimates

                # least squares fit of gaussian to data
                # p_hops = [np.mean(hop_lengths_boot), np.std(hop_lengths_boot)]
                # solp_hops, cov_x = curve_fit(gaussian, bins_hop_centered, hist_jump[i, :], p_hops,
                #                              bounds=([-np.inf, 0], [np.inf, np.inf]))
                #
                # mean_mle, sigma_mle = hop_mle
                # mean_ls, sigma_ls = solp_hops
                #
                # plt.bar(bins_hop_centered, hist_jump[i, :], hop_bin_width)
                # x = np.linspace(min(hop_lengths_boot), max(hop_lengths_boot), 500)
                # plt.plot(x, gaussian(x, mean_ls, sigma_ls), label='$\sigma$=%.3f nm (LS)' % sigma_ls, color='blue',
                #          linewidth=2)
                # plt.plot(x, gaussian(x, mean_mle, sigma_mle), label='$\sigma$=%.3f nm (MLE)' % sigma_mle, color='red',
                #          linewidth=2)
                # plt.legend(fontsize=14)
                # plt.xlabel('Hop length (z-direction, nm)', fontsize=14)
                # plt.ylabel('Probability', fontsize=14)
                # plt.gcf().get_axes()[0].tick_params(labelsize=14)
                # plt.tight_layout()
                # plt.show()
                # exit()

                self.hop_sigma_distribution.append(hop_mle[1])

                # Things for plotting
                if plot:

                    hop_mean += hop_mle[0]

                i += 1  # probably can go back to loop since using MLE
            except RuntimeError:  # sometimes bootstrapping gives unfittable data
                continue

        hop_mean /= nboot

        if plot:

            fig, ax = plt.subplots(2, 2, figsize=(12, 8))

            ax[0, 0].bar(bins_dwell_centered, hist_dwell.mean(axis=0), width=dwell_bin_width, align='center')
            ax[0, 0].set_xlabel('Dwell Time (ns)', fontsize=14)
            ax[0, 0].set_ylabel('Frequency', fontsize=14)

            ax[0, 0].plot(bins_dwell_centered, bins_dwell_centered ** -(1 + np.mean(self.alpha_distribution)),
                          '--', color='black', label=r'$\alpha_{fit}$ = %.3f $\pm$ %.3f ns$^{-1}$' %
                          (np.mean(self.alpha_distribution), np.std(self.alpha_distribution)))
            ax[0, 0].legend(fontsize=12)

            ax[0, 1].hist(self.alpha_distribution, bins=nbins)
            ax[0, 1].set_xlabel(r'Bootstrapped $\alpha$ (ns$^{-1}$)', fontsize=14)
            ax[0, 1].set_ylabel('Frequency', fontsize=14)

            ax[1, 0].bar(bins_hop_centered, hist_jump.mean(axis=0), width=hop_bin_width)
            ax[1, 0].plot(bins_hop_centered, gaussian(bins_hop_centered, hop_mean, np.mean(self.hop_sigma_distribution))
                          , '--', label='Gaussian fit\n $\sigma$=%.2f $\pm$ %.2f nm' %
                          (np.mean(self.hop_sigma_distribution), np.std(self.hop_sigma_distribution)), color='black')

            ax[1, 0].set_xlabel('Hop Length ($z$-direction, nm)', fontsize=14)
            ax[1, 0].set_ylabel('Frequency', fontsize=14)
            ax[1, 0].legend(fontsize=12)

            ax[1, 1].hist(self.hop_sigma_distribution, bins=nbins)
            ax[1, 1].set_xlabel('Bootstrapped $\sigma$ (nm)', fontsize=14)
            ax[1, 1].set_ylabel('Frequency', fontsize=14)

            plt.tight_layout()

            if save:
                plt.savefig('hop_dwell_distribution.pdf')

            if show:
                plt.show()

    def estimate_hurst(self, nboot=200, max_k=5):
        """ Estimate the hurst parameter by fitting the emperical autocovariance function to theory:

        \gamma(k) = \dfrac{1}{2}[|k-1|^{2H} - 2|k|^{2H} + |k + 1|^{2H}]

        :param nboot: number of bootstrap trials used to generate distribution of H's
        :param max_k: maximum number of time steps to fit in autocovariance function

        :type nboot: int
        :type max_k: int
        """

        max_hops = max([len(x) for x in self.hop_lengths])

        acf = np.zeros([len(self.hop_lengths), max_hops])

        keep = []  # list to hold indices of trajectories with a non-zero amount of hops
        for i in range(len(self.hop_lengths)):
            hops = self.hop_lengths[i]
            if len(hops) > 1:
                acf[i, :len(self.hop_lengths[i])] = timeseries.acf(self.hop_lengths[i])
                keep.append(i)

        acf = acf[keep, :]
        ntraj = len(keep)

        self.hop_acf = [acf[np.nonzero(acf[:, i]), i].mean() for i in range(max_hops)]

        plt.plot(self.hop_acf)
        plt.show()

        for b in range(nboot):

            traj = np.random.randint(0, ntraj, size=ntraj)

            hboot = acf[traj, 1].mean()

            H = np.log(2 * hboot + 2) / (2 * np.log(2))  # initial guess at H based on first dip in autocovariance

            # if H > 0:
            #     self.hurst_distribution.append(H)

            acf_boot = [acf[traj, i][np.nonzero(acf[traj, i])].mean() for i in range(max_k + 1)]

            h_opt = curve_fit(fitting_functions.hurst_autocovariance, np.arange(max_k + 1), acf_boot[:(max_k + 1)],
                              p0=H)[0]

            if h_opt > 0:

                self.hurst_distribution.append(h_opt[0])

            #self.hurst_distribution.append(np.log(2 * hboot + 2) / (2 * np.log(2)))

    def update_database(self, file="msd.db", tablename="msd", type='parameters', data=None):
        """ Update SQL database with information from this run

        :param file: relative path (relative to directory where this script is stored) to database to be updated
        :param tablename: name of table being modified in database
        :param type: The type of info to be updated/added to the table. 'parameters' indicates an update to alpha,
        sigma, hurst, sim_length and mw. 'msds' indicates an update to python_MSD, python_MSD_CI_upper and
        python_MSD_CI_Lower
        :param data: data to be filled in if it is not a part of this object. For example, to update python MSDs,
        include a list with [python_MSD, python_MSD_CI_lower, python_MSD_CI_upper] in that order.

        :type file: str
        :type tablename: str
        :type type: str
        :type data: list
        """

        connection = sql.connect("%s/%s" % (self.location, file))
        crsr = connection.cursor()

        alpha = np.mean(self.alpha_distribution)
        sigma = np.mean(self.hop_sigma_distribution)
        hurst = np.mean(self.hurst_distribution)

        check_existence = "SELECT COUNT(1) FROM %s WHERE name = '%s' and penalty = %.2f" % (tablename,
                           self.residue.name, self.breakpoint_penalty)

        output = crsr.execute(check_existence).fetchall()

        if type == 'parameters':

            if output[0][0] == 1:

                update_entry = "UPDATE %s SET alpha = %.2f, sigma = %.2f, hurst = %.2f, sim_length = %.2f, mw = %.2f " \
                               "WHERE name = '%s' and penalty = %.2f" % (tablename, alpha, sigma, hurst, self.time[-1],
                                                                         self.residue.mw, self.residue.name,
                                                                         self.breakpoint_penalty)

                crsr.execute(update_entry)

            else:

                fill_new_entry = "INSERT INTO %s (name, alpha, sigma, hurst, penalty, sim_length, mw) VALUES ('%s', " \
                                 "%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)" % (tablename, self.residue.name, alpha, sigma,
                                                                          hurst, self.breakpoint_penalty, self.time[-1],
                                                                          self.residue.mw)

                crsr.execute(fill_new_entry)

        elif type == 'msd_ensemble' or type == 'msd_time_average':

            if type == 'msd_ensemble':
                data_labels = ['python_MSD', 'python_MSD_CI_lower', 'python_MSD_CI_upper']
            elif type == 'msd_time_average':
                data_labels = ['python_TAMSD', 'python_TAMSD_CI_lower', 'python_TAMSD_CI_upper']

            error_message = "You must provide the MSD data to be input into the database in the format " \
                                    "[python_MSD, python_MSD_CI_lower, python_MSD_CI_upper]"

            if output[0][0] > 0:

                try:
                    update_entry = "UPDATE %s SET %s = %.3f, %s = %.3f, %s = %.3f WHERE name = '%s' and penalty = %.2f" \
                                   % (tablename, data_labels[0], data[0], data_labels[1], data[1], data_labels[2],
                                      data[2], self.residue.name, self.breakpoint_penalty)
                except TypeError:
                    raise TypeError(error_message)
                except IndexError:
                    raise IndexError(error_message)

                crsr.execute(update_entry)

            else:

                try:
                    fill_new_entry = "INSERT INTO %s (name, %s, %s, %s, penalty) VALUES ('%s', %.3f, %.3f, %.3f, %.2f)" \
                                     % (tablename, data_labels[0], data_labels[1], data_labels[2], self.residue.name,
                                     data[0], data[1], data[2], self.breakpoint_penalty)
                except TypeError:
                    raise TypeError(error_message)
                except IndexError:
                    raise IndexError(error_message)

                crsr.execute(fill_new_entry)
        else:
            sys.exit('Please choose a ')

        connection.commit()
        connection.close()


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.load:

        sys = file_rw.load_object(args.load)

    else:

        sys = System(args.trajectory, args.gro, args.residue, start=args.begin, end=args.end, step=args.step,
                     ma=args.moving_average)

        sys.calculate_solute_partition(spline=False, membrane_residue='HII')

        sys.hops_and_dwells(penalty=args.breakpoint_penalty)

        sys.fit_distributions(nbins=args.nbins, nboot=args.nboot, plot=True, show=True)

        sys.estimate_hurst()

        if args.update:
            sys.update_database()

        file_rw.save_object(sys, 'forecast_%s.pl' % args.residue)

    sys.estimate_hurst()
    plt.hist(sys.hurst_distribution, bins=25)
    plt.show()
    exit()
    # plt.plot(sys.hop_acf, linewidth=2,  label='Simulation')
    # plt.xlabel('k', fontsize=14)
    # plt.ylabel('Autocovariance', fontsize=14)
    # plt.plot(fitting_functions.hurst_autocovariance(np.arange(20), np.mean(sys.hurst_distribution)), label='Analytical', linewidth=2)
    # plt.gcf().get_axes()[0].tick_params(labelsize=14)
    # plt.legend(fontsize=14)
    # plt.xlim(-1, 20)
    # plt.show()

    # sys.location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # This script location
    # sys.breakpoint_penalty = 0.25

    # sys.estimate_hurst()
    # sys.hops_and_dwells(penalty=1)
    # file_rw.save_object(sys, 'forecast_%s.pl' % args.residue)

    # histogram of time spent in pore region. Kind of interesting
    # plt.hist(np.sum(sys.partition, axis=0) / sys.partition.shape[0], bins=50)
    # plt.show()

    # simulate ntrajsim trajectories for same length as MD
    random_walks = CTRW(2000, args.ntrajsim, dt=sys.dt, hop_dist='fbm', dwell_dist='power')
    # sys.hurst_distribution = 0.5*np.ones(10)
    # sys.alpha_distribution = np.ones(10)
    random_walks.generate_trajectories(fixed_time=True, distributions=(sys.alpha_distribution,
                                       sys.hop_sigma_distribution, sys.hurst_distribution), discrete=True, ll=1)
    random_walks.calculate_msd(ensemble=args.ensemble)

    if args.ensemble:  # Ensemble-averaged MSD

        random_walks.bootstrap_msd(fit_power_law=True)
        random_walks.plot_msd(plot_power_law=True, show=False)
        # sys.update_database(type='msd', data=[1000 * (x / sys.time[-1]) for x in random_walks.final_msd])
        if args.update:
            sys.update_database(type='msd_ensemble', data=random_walks.final_msd)

    else:  # Time-averaged MSD

        random_walks.bootstrap_msd(fit_linear=False)
        random_walks.plot_msd(show=False, end_frame=8000)  # get data up to 400 ns
        if args.update:
            sys.update_database(type='msd_time_average', data=random_walks.final_msd)
