#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import ruptures
from scipy.optimize import curve_fit, minimize
from LLC_Membranes.llclib import topology, physical, file_rw, fitting_functions, timeseries
import tqdm
import sqlite3 as sql
import os
import sys


def gaussian(points, mean, sigma):

    return (1 / np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(
        -(points - mean) ** 2 / (2 * sigma ** 2))


class SFBMParameters(object):

    def __init__(self, traj, gro, res, start=0, end=-1, step=1, ma=False, nmodes=1):
        """ Using an MD trajectory, calculate fit parameters for the dwell time distribution and hop length distribution

        :param traj: name of MD trajectory (GROMACS .xtc or .trr)
        :param gro: name of .gro (or .pdb) coordinate file associated with same topology as traj
        :param res: name of residue that is hopping around
        :param start: first frame of trajectory to include
        :param end: last frame of trajectory to include
        :param step: include every step frames in calculations
        :param ma: calculate a moving average on the center of mass coordinate traces
        :param nmodes: number of modes to consider. 1 mode will treat all hops and dwells as being performed in the
        same environment. 2 modes will treat hops and dwells in the tails versus pores differently.

        :type traj: str
        :type gro: str
        :type res: str
        :type start: int
        :type end: int
        :type step: int
        :type nmodes: int
        """

        # load trajectory
        print('Loading trajectory...', end='', flush=True)
        self.t = md.load(traj, top=gro)[start:end:step]
        print('Done!')

        self.time = self.t.time / 1000  # time in nanoseconds
        self.dt = self.time[1] - self.time[0]  # time step
        self.nT = self.t.n_frames

        keep = [a.index for a in self.t.topology.atoms if a.residue.name == res]  # get indices of atoms making up residue of interest

        self.res_start = keep[0]  # index where residue starts

        self.residue = topology.Residue(res)
        self.nres = len(keep) // self.residue.natoms  # number of residues in system
        self.mass = [v for v in self.residue.mass.values()]  # mass of atoms making up residue
        self.pos = self.t.xyz[:, keep, :]  # positions of residue atoms

        print('Calculating centers of mass...', end='', flush=True)
        self.com = physical.center_of_mass(self.pos, self.mass)  # center of mass of residues
        print('Done!')

        # plot z-coordinate trace
        # for i in range(self.com.shape[1]):
        #     plt.plot(self.time, self.com[:, i, 2], linewidth=2)
        #     plt.tick_params(labelsize=14)
        #     plt.xlabel('Time (ns)', fontsize=14)
        #     plt.ylabel('$z$-coordinate (nm)', fontsize=14)
        #     plt.tight_layout()
        #     plt.show()

        # plot first order difference histogram
        # nbins = 100
        # plt.hist((self.com[1:, :, 2] - self.com[:-1, :, 2]).flatten(), bins=nbins, density=True)
        # plt.tick_params(labelsize=14)
        # plt.xlabel('$z$-direction hop length (nm)', fontsize=14)
        # plt.ylabel('Frequency', fontsize=14)
        # plt.tight_layout()
        # plt.show()
        # exit()

        #plot first order differences
        # plt.figure()
        # plt.plot((self.com[1:, 1, 2] - self.com[:-1, 1, 2]))
        # # plt.plot(self.com[:, 1, 2], linewidth=2)
        # plt.tick_params(labelsize=14)
        # plt.xticks(np.linspace(0, 2000, 9), [int(i) for i in np.linspace(0, 1000, 9)])
        # plt.xlabel('Time (ns)', fontsize=14)
        # plt.ylabel('$z$-direction hop length (nm)', fontsize=14)
        # plt.tight_layout()
        # plt.show()
        # exit()

        #create timeseries by randomly drawing from unconditional pdf of hop lengths
        # plt.figure()
        # p = (self.com[1:, :, 2] - self.com[:-1, :, 2]).flatten()
        # for i in range(10):
        #     t = np.random.choice(p, size=2000)
        #     plt.plot(np.linspace(0, 1000, len(t)), np.cumsum(t), linewidth=2)
        #     plt.tick_params(labelsize=14)
        #     plt.xlabel('Time (ns)', fontsize=14)
        #     plt.ylabel('$z$-coordinate (nm)', fontsize=14)
        #     plt.tight_layout()
        #     plt.show()
        # exit()

        if ma:
            self.calculate_moving_average(ma)

        # initialize for later
        self.nmodes = nmodes
        self.pore_centers = None
        self.dwell_times = [[] for _ in range(self.nmodes)]  # dwell times in each dynamical mode
        self.tail_dwells = []  # tail as in tail ends of trajectory
        self.hop_lengths = []
        self.hurst_distribution = [[] for _ in range(self.nmodes)]
        self.hop_acf = None
        self.alpha_distribution = [[] for _ in range(self.nmodes)]  # distribution of alpha for poisson process
        self.hop_sigma_distribution = [[] for _ in range(self.nmodes)]  # distribution of standard deviation of hop lengths
        self.breakpoint_penalty = 0
        self.location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # This script location
        self.r = []  # radial distance of solute from pore center when hops are made
        self.segments = [[] for _ in range(self.nmodes)]
        self.hop_series = [[] for _ in range(self.nmodes)]
        self.nsolute = self.com.shape[1]
        self.transition_matrix = np.zeros([self.nmodes, self.nmodes])

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
        :param write_tcl: write a tcl script that will color code solutes based on their radial position. A good way to
        make sure this function worked as intended.

        :type r: float
        :type spline: bool
        :type membrane_residue: str
        :type write_tcl: bool
        """

        membrane = topology.LC('%s' % membrane_residue)  # object w/ attributes of LC making up membrane
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

        :param frame: frame to write .tcl script for. Make sure you have the .gro file corresponding to this frame or
        the output will not be useful. (use get_frame.py if you don't have the proper gro file)
        :param name: name of .tcl file

        :type frame: int
        :type name: str
        """

        pore_indices = [self.res_start + self.residue.natoms * i for i in np.where(self.partition[frame, :])[0]]
        tail_indices = [self.res_start + self.residue.natoms * i for i in np.where(self.partition[frame, :] == False)[0]]  # have to use double equals sign. Using is doesn't work with np.where

        with open(name, 'w') as f:
            f.write('color Display Background white\n')
            f.write('mol addrep 0\n')
            f.write('mol modselect 0 0 index')
            for i in pore_indices:
                end = i + len(self.mass)
                f.write(' %s to %s' % (i, end - 1))
            f.write('\n')
            f.write('mol modcolor 0 0 ColorID 0\n')
            f.write('mol modstyle 0 0 CPK 2.0 0.3 12.0 12.0\n')
            f.write('mol addrep 0\n')
            f.write('mol modselect 1 0 index')
            for i in tail_indices:
                end = i + len(self.mass)
                f.write(' %s to %s' % (i, end - 1))
            f.write('\n')
            f.write('mol modstyle 1 0 CPK 2.0 0.3 12.0 12.0\n')
            f.write('mol modcolor 1 0 ColorID 1\n')

    def write_spline_coordinates(self, frame=-1):
        """ write spline coordinates to a .gro file
        """

        coordinates = np.reshape(self.pore_centers[frame, ...], (self.pore_centers.shape[1]*self.pore_centers.shape[2],
                                                                 3))
        file_rw.write_gro_pos(coordinates, 'spline.gro', name='K')

    def hops_and_dwells(self, penalty=0.25, locations=False):
        """ Find breakpoints then assemble lists of dwell times and hop lengths. See documentation for Ruptures python
        package: http://ctruong.perso.math.cnrs.fr/ruptures-docs/build/html/index.html for more options that can be
        added

        :param penalty: penalty for cost function minimization
        :param nframes_dwell: number of frames that a solute needs to stay in the region of interest in order to be
        analyzed for hops.
        :param locations: if True, this function will return the distance from the pore center when hops are made

        :type penalty: float
        :type nframes_dwell: int
        :type locations: bool
        """

        self.breakpoint_penalty = penalty

        # Handle initial and final dwell times since they don't have a true beginning or end. If the dwell times are
        # longer than some limit, then they will be included so that we don't miss out on dwells that last on the
        # order of the entire trajectory length.
        initial_dwells = [[] for _ in range(self.nmodes)]  # length time particle initially dwells before its first hop
        final_dwells = [[] for _ in range(self.nmodes)]  # length of time particle dwells up until end of trajectory

        for j in tqdm.tqdm(range(self.com.shape[1])):

            hop_sequence = []  # create sequence of hop lengths IF the mode does not change from previous mode

            hop_lengths = [[] for _ in range(self.nmodes)]

            # find all hops
            bp = ruptures.detection.Binseg(model='l2', min_size=1, jump=1).fit_predict(self.com[:, j, :]
                                                                                       , pen=self.breakpoint_penalty)

            # ruptures.display(self.com[:, j, 2], bp)
            # plt.xlabel('Time (ns)', fontsize=14)
            # plt.ylabel('$z$-position (nm)', fontsize=14)
            # plt.xticks(np.linspace(0, 2000, 9), [int(i) for i in np.linspace(0, 1000, 9)])
            # plt.tight_layout()
            # plt.show()
            # exit()

            if locations:  # hasn't been tested with modified script

                pore = j % self.pore_centers.shape[1]  # assumes solutes placed with place_solutes_pores.py

                if len(self.pore_centers.shape) == 4:  # spline

                    pore_centers = self.pore_centers[:, pore, ...][bp[:-1]]
                    coms = self.com[:, j, :][bp[:-1]]
                    box = self.t.unitcell_vectors[:, ...][bp[:-1]]

                    for k in range(len(bp[:-1])):
                        self.r += physical.radial_distance_spline(pore_centers[k, ...], coms[k, :][np.newaxis, :],
                                                                  box[k]).tolist()

                else:  # no spline

                    pore_centers = self.pore_centers[:, pore, :][bp[:-1]]
                    coms = self.com[:, j, :2][bp[:-1]]
                    self.r += np.linalg.norm(pore_centers - coms, axis=1).tolist()

            location = 0
            if self.nmodes == 2:
                location = self.partition[bp[0], j]

            initial_dwells[location].append(bp[0])
            if len(bp) > 1:  # handle case where the only break point is at the end of the simulation
                if self.nmodes == 2:
                    location = self.partition[bp[-2], j]  # -1 is always the last frame
                final_dwells[location].append(bp[-1] - bp[-2])

            previous_location = -1
            for k in range(len(bp) - 1):

                if k > 0:  # need two different z coordinates to get a hop length
                    hop_length = np.mean(self.com[bp[k]:bp[k + 1], j, 2]) - np.mean(self.com[bp[k - 1]:bp[k], j, 2])
                else:
                    # bp[0] != 0, rather it is the first break point location. bp[-1] is the last point in the
                    # time series though
                    # hop at first segment
                    hop_length = np.mean(self.com[bp[k]:bp[k + 1], j, 2]) - np.mean(self.com[0:bp[k], j, 2])

                location = 0
                if self.nmodes == 2:
                    location = self.partition[bp[k], j]  # determine location of solute when it hops

                hop_lengths[location].append(hop_length)

                if location == previous_location:
                    hop_sequence.append(hop_length)
                else:
                    if hop_sequence:  # avoid appending empty lists. This is possible since I reset the list each traj
                        self.hop_series[previous_location].append(hop_sequence)
                    hop_sequence = [hop_length]

                # Assume solute stays in region it hopped to
                if (k + 1) < (len(bp) - 1):  # exclude first and last segments for dwell times
                    self.dwell_times[location].append(bp[k + 1] - bp[k])  # discrete form

                previous_location = location

            # make sure to get last observed sequence
            if hop_sequence:
                self.hop_series[previous_location].append(hop_sequence)
            #
            # # Analyze sub-trajectories where solute is in defined pore region
            # begin = 0
            # for i, end in enumerate(switch_points[::2]):  # only at switch points after which solute is in pores
            #     if (end - begin) > nframes_dwell:  # only include hops/dwells from trajectories long enough to analyze
            #         # np.savez_compressed('MET_trace.npz', z=self.com[begin:end, j, :])
            #         # exit()
            #         # movement_3d = np.linalg.norm(self.com[begin:end, j, :], axis=1)  # magnitude of distance travelled
            #         bp = ruptures.detection.Binseg(model='l2', min_size=1, jump=1).fit_predict(self.com[begin:end, j, :]
            #                                        , pen=self.breakpoint_penalty)
            #         # if j == 1:
            #         #     # bp = ruptures.detection.Binseg(model='l2', min_size=1, jump=1).fit_predict(self.com[begin:end, j, 2], pen=self.breakpoint_penalty)
            #         #     ruptures.display(self.com[begin:end, j, 2], bp)
            #         #     plt.xlabel('Time (ns)', fontsize=14)
            #         #     plt.ylabel('$z$-position (nm)', fontsize=14)
            #         #     plt.xticks(np.linspace(0, 2000, 9), [int(i) for i in np.linspace(0, 1000, 9)])
            #         #     plt.tight_layout()
            #         #     plt.show()
            #         #     exit()

                    # for visualization predicted z-hops on top of time series
                    # traj_hops = np.zeros([2*len(bp) - 1, 2])
                    #
                    # traj_hops[1::2, 0] = self.time[bp[:-1]]*2
                    # traj_hops[2::2, 0] = self.time[bp[:-1]]*2

                        # traj_hops[2*k, 0] = bp[k]
                        # traj_hops[2*k + 1, 0] = bp[k + 1]
                        # traj_hops[2*k:(2*k + 2), 1] = np.mean(self.com[bp[k]:bp[k + 1], j, 2])

                    # ruptures.display(self.com[begin:end, j, :], bp, figsize=(10, 6))
                    # plt.plot(traj_hops[:, 0], traj_hops[:, 1], '--', color='black')
                    # plt.show()
                    # exit()
                # try:
                #     begin = switch_points[2*i + 1]
                # except IndexError:
                #     pass

            self.hop_lengths.append(hop_lengths)

        # tail_dwells as in tail ends of trajectory
        self.tail_dwells = [initial_dwells[i] + final_dwells[i] for i in range(self.nmodes)]

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
        hist_dwell = np.zeros([self.nmodes, nboot, nbins])
        hist_jump = np.zeros([self.nmodes, nboot, nbins])

        all_hops = [[] for _ in range(self.nmodes)]

        for h in self.hop_lengths:
            for m in range(self.nmodes):
                all_hops[m] += h[m]

        hop_mean = np.zeros(self.nmodes)

        # Define bins so histograms can be added together. Do it separately for each mode
        min_dwells = [min(i) for i in self.dwell_times]
        max_dwells = [max(max(self.dwell_times[i]), max(self.tail_dwells[i])) for i in range(self.nmodes)]
        min_hops = [min(i) for i in all_hops]
        max_hops = [max(i) for i in all_hops]

        bins_dwell = np.zeros([self.nmodes, nbins + 1])
        bins_hop = np.zeros([self.nmodes, nbins + 1])
        for m in range(self.nmodes):
            bins_dwell[m, :] = np.linspace(min_dwells[m], max_dwells[m], nbins + 1)
            bins_hop[m, :] = np.linspace(min_hops[m], max_hops[m], nbins + 1)

        dwell_bin_width = bins_dwell[:, 1] - bins_dwell[:, 0]
        hop_bin_width = bins_hop[:, 1] - bins_hop[:, 0]

        bins_dwell_centered = np.zeros([self.nmodes, nbins])
        bins_hop_centered = np.zeros([self.nmodes, nbins])
        for m in range(self.nmodes):
            bins_dwell_centered[m, :] = [i + dwell_bin_width[m]/2 for i in bins_dwell[m, :-1]]
            bins_hop_centered[m, :] = [i + hop_bin_width[m]/2 for i in bins_hop[m, :-1]]

        i = 0
        while i < nboot:

            print('Bootstrapping ... Trial %s\r' % (i + 1), end='', flush=True)
            # dwell times
            for m in range(self.nmodes):

                try:

                    # Recreate dwell time distribution by randomly choosing from all dwell times with replacement
                    dwell_times_boot = np.random.choice(self.dwell_times[m], size=len(self.dwell_times[m]),
                                                        replace=True)

                    # Add long dwell times randomly selected from beginning and end of trajectory
                    tail_dwells = np.array(self.tail_dwells[m])[np.where(self.tail_dwells[m] > max(dwell_times_boot))[0]]

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

                    self.alpha_distribution[m].append(maximum - 1)

                    # Graphical representation of log-likelihood maximization
                    # ll = []  # list of log-likelihoods for each alpha tested
                    # alphas = np.linspace(1.01, 3, 100)
                    # for a in alphas:
                    #     ll.append(fitting_functions.power_law_discrete_log_likelihood(a, dwell_times_boot,
                    #                                                                   int(min(dwell_times_boot))))
                    #
                    # plt.plot(alphas, ll)
                    # plt.show()
                    # exit()

                    if plot:

                        # bin the bins
                        hist_dwell[m, i, :], _ = np.histogram(dwell_times_boot, bins_dwell[m, :], density=True)

                    # jump lengths
                    hop_lengths_boot = np.random.choice(all_hops[m], size=len(all_hops[m]), replace=True)

                    if plot:

                        hist_jump[m, i, :], _ = np.histogram(hop_lengths_boot, bins_hop[m, :], density=True)

                    # Maximum likelihood estimate of mean and sigma of hop length distribution
                    hop_mle = [np.mean(hop_lengths_boot), np.std(hop_lengths_boot)]  # mean and std are mle estimates

                    self.hop_sigma_distribution[m].append(hop_mle[1])

                    # Things for plotting
                    if plot:

                        hop_mean += hop_mle[0]

                    if m == self.nmodes - 1:
                        i += 1  # probably can go back to loop since using MLE

                except RuntimeError:  # sometimes bootstrapping gives unfittable data
                    continue

        hop_mean /= nboot

        if plot:

            fig, ax = plt.subplots(self.nmodes, 2, figsize=(self.nmodes * 6, 8))

            if self.nmodes == 1:
                ax = ax[np.newaxis, :]
            # hops, hop_ax = plt.subplots()
            # dwells, dwell_ax = plt.subplots()
            fontsize = 16
            label_fontsize = 14
            legend_fontsize = 12

            for m in range(self.nmodes):

                ax[m, 0].bar(bins_dwell_centered[m, :], hist_dwell[m, ...].mean(axis=0), width=dwell_bin_width[m],
                             align='center')
                ax[m, 0].set_xlabel('Dwell Time (ns)', fontsize=label_fontsize)
                ax[m, 0].set_ylabel('Frequency', fontsize=label_fontsize)

                ax[m, 0].plot(bins_dwell_centered[m, :], bins_dwell_centered[m, :] ** -(1 +
                              np.mean(self.alpha_distribution)), '--', color='black',
                              label=r'$\alpha_{fit}$ = %.3f $\pm$ %.3f ns$^{-1}$' % (np.mean(self.alpha_distribution[m]),
                                                                                     np.std(self.alpha_distribution[m])))
                ax[m, 0].legend(fontsize=legend_fontsize)
                ax[m, 0].tick_params(labelsize=fontsize)

                ax[m, 1].bar(bins_hop_centered[m, :], hist_jump[m, ...].mean(axis=0), width=hop_bin_width[m])
                ax[m, 1].plot(bins_hop_centered[m, :],
                              gaussian(bins_hop_centered[m, :], hop_mean[m], np.mean(self.hop_sigma_distribution[m]))
                              , '--', label='Gaussian fit\n $\sigma$=%.2f $\pm$ %.2f nm' %
                              (np.mean(self.hop_sigma_distribution[m]), np.std(self.hop_sigma_distribution[m])),
                              color='black')
                ax[m, 1].legend(fontsize=legend_fontsize)
                ax[m, 1].set_xlabel('Hop Length (nm)', fontsize=label_fontsize)
                ax[m, 1].set_ylabel('Frequency', fontsize=label_fontsize)

                ax[m, 1].tick_params(labelsize=fontsize)

            # dwell_ax.bar(bins_dwell_centered, hist_dwell.mean(axis=0), width=dwell_bin_width, align='center')
            # dwell_ax.set_xlabel('Dwell Time (ns)', fontsize=fontsize)
            # dwell_ax.set_ylabel('Frequency', fontsize=fontsize)
            #
            # dwell_ax.plot(bins_dwell_centered, bins_dwell_centered ** -(1 + np.mean(self.alpha_distribution)),
            #               '--', color='black', label='Power law MLE fit', linewidth=2)
            # dwell_ax.legend(fontsize=fontsize)
            # dwell_ax.text(290, 0.022, r'$p \propto t^{(-1 - \alpha)}$', fontsize=fontsize)
            # dwell_ax.text(290, 0.019, r'$\alpha_{fit}$ = %.3f $\pm$ %.3f ns$^{-1}$' %
            # (np.mean(self.alpha_distribution), np.std(self.alpha_distribution)), fontsize=fontsize)
            # dwell_ax.tick_params(labelsize=fontsize)

            # ax[0, 1].hist(self.alpha_distribution, bins=nbins)
            # ax[0, 1].set_xlabel(r'Bootstrapped $\alpha$ (ns$^{-1}$)', fontsize=14)
            # ax[0, 1].set_ylabel('Frequency', fontsize=14)
            #
            # ax[1, 0].bar(bins_hop_centered, hist_jump.mean(axis=0), width=hop_bin_width)
            # ax[1, 0].plot(bins_hop_centered, gaussian(bins_hop_centered, hop_mean, np.mean(self.hop_sigma_distribution))
            #               , '--', label='Gaussian fit\n $\sigma$=%.2f $\pm$ %.2f nm' %
            #               (np.mean(self.hop_sigma_distribution), np.std(self.hop_sigma_distribution)), color='black')
            #
            # ax[1, 0].set_xlabel('Hop Length ($z$-direction, nm)', fontsize=14)
            # ax[1, 0].set_ylabel('Frequency', fontsize=14)
            # ax[1, 0].legend(fontsize=12)
            #
            # hop_ax.bar(bins_hop_centered, hist_jump.mean(axis=0), width=hop_bin_width)
            # hop_ax.plot(bins_hop_centered, gaussian(bins_hop_centered, hop_mean, np.mean(self.hop_sigma_distribution))
            #               , '--', label='Gaussian MLE fit', color='black', linewidth=2)
            #
            # hop_ax.set_xlabel('Hop Length ($z$-direction, nm)', fontsize=fontsize)
            # hop_ax.set_ylabel('Frequency', fontsize=fontsize)
            # hop_ax.legend(fontsize=fontsize)
            # hop_ax.text(0.31, 0.95, '$p \propto e^{\dfrac{(x-\mu)^2}{2\sigma^2}}$', fontsize=fontsize)
            # hop_ax.text(0.31, 0.825, '$\sigma$=%.2f $\pm$ %.2f nm' % (np.mean(self.hop_sigma_distribution),
            #                                                     np.std(self.hop_sigma_distribution)), fontsize=fontsize)
            # hop_ax.tick_params(labelsize=fontsize)
            #
            # ax[1, 1].hist(self.hop_sigma_distribution, bins=nbins)
            # ax[1, 1].set_xlabel('Bootstrapped $\sigma$ (nm)', fontsize=14)
            # ax[1, 1].set_ylabel('Frequency', fontsize=14)

            fig.tight_layout()
            # hops.tight_layout()
            # dwells.tight_layout()

            if save:
                fig.savefig('hop_dwell_distribution.pdf')
                # hops.savefig('hop_distribution.pdf')
                # dwells.savefig('dwell_distribution.pdf')

            if show:
                plt.show()

    def identify_segments(self, min_segment_length=10):
        """ Identify all segments between region transitions that are at least min_segment_length long
        TODO: add this to identify_states.py -- can be used to calculate acf in each region
        :param min_segment_length: minimum number of frames for a segment to be kept

        :type min_segment_length: int

        :return max_dwell: the maximum segment length
        """

        nsolute = self.com.shape[1]  # number of solute trajectories
        max_dwell = 0

        for s in range(nsolute):

            # Get indices where res jumps between regions
            switch_points = timeseries.switch_points(self.partition[:, s])
            # switch_points = np.argwhere(np.diff(self.partition[:, s])).squeeze().tolist()
            #
            # # add last frame as a switch point
            # try:
            #     switch_points.append(self.partition.shape[0])
            # except AttributeError:  # if there are no switches, it won't return a list
            #     switch_points = list(switch_points)
            #     switch_points.append(self.partition.shape[0])
            #
            # switch_points = np.array([0] + switch_points)  # also add first frame

            dwells = switch_points[1:] - switch_points[:-1]  # length of dwells

            if max(dwells) > max_dwell:
                max_dwell = max(dwells)

            long_enough = np.where(dwells >= 10)[0]  # indices of segments with dwell times that are long enough

            for seg in long_enough:  # make pairs of indices (beginning and end) defining where eligible segments are
                segment = (switch_points[seg], switch_points[seg + 1])
                location = self.partition[segment[0]:segment[1], s]
                if self.nmodes > 1:  # if this is a multi-mode model, identify the mode and assign segment accordingly
                    self.segments[int(round(location.mean()))].append(segment)
                else:
                    self.segments.append(segment)

        return max_dwell

    def estimate_hurst(self, nboot=200, max_k=5, confidence=95):
        """ Estimate the hurst parameter by fitting the emperical autocovariance function to theory:

        \gamma(k) = \dfrac{1}{2}[|k-1|^{2H} - 2|k|^{2H} + |k + 1|^{2H}]

        :param nboot: number of bootstrap trials used to generate distribution of H's
        :param max_k: maximum number of time steps to fit in autocovariance function
        :param confidence: confidence interval

        :type nboot: int
        :type max_k: int
        :type confidence: float
        """

        # lists of arrays that will be needed
        self.hop_acf = []  # average autocorrelation function
        lower_confidence = (100 - confidence) / 2  # confidence intervals (percentiles)
        upper_confidence = 100 - lower_confidence

        # plot settings. Work around for plot subscripting
        if self.nmodes == 2:
            fig, (ax1, ax2) = plt.subplots(1, self.nmodes, figsize=(6 * self.nmodes, 5))
            ax = [ax1, ax2]
        else:
            fig, ax1 = plt.subplots(1, 1, figsize=(6, 5))
            ax = [ax1]

        fontsize = 14

        for m in range(self.nmodes):

            nhops = np.array([len(x) for x in self.hop_series[m]])

            max_hops = max(nhops)

            acf = np.zeros([len(self.hop_series[m]), max_hops])
            boot = np.zeros([nboot, max_hops])

            keep = []  # list to hold indices of trajectories with a non-zero amount of hops
            for i in range(acf.shape[0]):
                hops = self.hop_series[m][i]
                if len(hops) > 2:  # correlation between two points is useless. Will always be +1, -1
                    autocorrelation = timeseries.acf(hops)
                    acf[i, :autocorrelation.size] = autocorrelation
                    keep.append(i)

            acf = acf[keep, :]
            nhops = nhops[keep]

            nsegments = acf.shape[0]

            for b in range(nboot):
                sol = np.random.randint(acf.shape[0], size=acf.shape[0])
                for i in range(max(nhops[sol])):
                    ndx = sol[np.nonzero(acf[sol, i])]
                    if not list(ndx):
                        boot[b, i] = 0
                    else:
                        boot[b, i] = acf[ndx, i].mean()
                    # This is more pythonic, but multiple exceptions are raised so it doesn't really work as intended.
                    # try:
                    #     boot[b, i] = acf[ndx, i].mean()
                    # except RuntimeWarning:  # happens if the solute with the max_hops dwell time is not included in 'sol'
                    #     boot[b, i] = 0

            # boot[b, :] = [acf[np.nonzero(acf[sol, i]), i].mean() for i in range(max_hops)]
            # bootstrap acf
            #self.hop_acf.append([acf[np.nonzero(acf[:, i]), i].mean() for i in range(max_hops)])
            self.hop_acf.append(boot.mean(axis=0))

            errorbars = np.zeros([2, len(self.hop_acf[m])])
            errorbars[0, :] = np.abs(np.percentile(boot, lower_confidence, axis=0) -
                                     boot.mean(axis=0))  # 2.5 percent of data below this value
            errorbars[1, :] = np.percentile(boot, upper_confidence, axis=0) - boot.mean(axis=0)

            for b in range(nboot):

                traj = np.random.randint(0, nsegments, size=nsegments)

                hboot = acf[traj, 1].mean()  # first time lag autocovariance

                H = np.log(2 * hboot + 2) / (2 * np.log(2))  # initial guess at H based on first dip in autocovariance

                # if H > 0:
                #     self.hurst_distribution.append(H)
                # can probably combine this loop with above bootstrapping loop
                acf_boot = [acf[traj, i][np.nonzero(acf[traj, i])].mean() for i in range(max_k + 1)]

                h_opt = curve_fit(fitting_functions.hurst_autocovariance, np.arange(max_k + 1), acf_boot[:(max_k + 1)],
                                  p0=H)[0]

                if h_opt > 0:

                    self.hurst_distribution[m].append(h_opt[0])

            ax[m].plot(self.hop_acf[m], linewidth=2, label='Simulated autocorrelation')
            ax[m].plot(np.arange(max_hops), fitting_functions.hurst_autocovariance(np.arange(max_hops),
                     np.mean(self.hurst_distribution[m])), '--', color='black', linewidth=2,
                     label='Fit theoretical autocorrelation')
            ax[m].fill_between(np.arange(max_hops), errorbars[1, :] + boot.mean(axis=0),
                               boot.mean(axis=0) - errorbars[0, :], alpha=0.25)
            ax[m].text(3.0, 0.5, '$\gamma(k) = \dfrac{1}{2}[|k-1|^{2H} - 2|k|^{2H} + |k + 1|^{2H}]$', fontsize=fontsize)
            ax[m].text(3.0, 0.3, '$H$ = %.2f $\pm$ %.2f' % (np.mean(self.hurst_distribution[m]),
                                                            np.std(self.hurst_distribution[m])), fontsize=fontsize)
            ax[m].set_xticks([1, 3, 5, 7, 9, 11, 13, 15])
            ax[m].set_xlabel('Lag (k)', fontsize=fontsize)
            ax[m].set_ylabel('Autocorrelation', fontsize=fontsize)
            ax[m].tick_params(labelsize=fontsize)
            ax[m].set_xlim(0.05, 15)
            ax[m].set_ylim(-1, 1)
            ax[m].legend(fontsize=fontsize)

        fig.tight_layout()
        fig.savefig('hop_acf.pdf')
        plt.show()

    def determine_transition_matrix(self, start=1):
        """ Create a transition matrix describing the probability of transitions between states. This is same as
        make_transition_matrix in identify_states.py. Could be added to a library. Would need to recast self.partition
        as an integer array.

        :param start: first frame to include. Must be greater than 0 since we will need to know the previous state

        :type start: int
        """

        nT = self.partition.shape[0]

        count_matrix = np.zeros_like(self.transition_matrix, dtype=int)

        for t in range(start, nT):  # start at frame 1. May need to truncate more as equilibration
            transitioned_from = self.partition[t - 1, :].astype(int)
            transitioned_to = self.partition[t, :].astype(int)
            for pair in zip(transitioned_from, transitioned_to):
                count_matrix[pair[0], pair[1]] += 1

        # normalize so rows sum to unity
        self.transition_matrix = (count_matrix.T / count_matrix.sum(axis=1)).T

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
                                                                         self.residue.MW, self.residue.name,
                                                                         self.breakpoint_penalty)

                crsr.execute(update_entry)

            else:

                fill_new_entry = "INSERT INTO %s (name, alpha, sigma, hurst, penalty, sim_length, mw) VALUES ('%s', " \
                                 "%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)" % (tablename, self.residue.name, alpha, sigma,
                                                                          hurst, self.breakpoint_penalty, self.time[-1],
                                                                          self.residue.MW)

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
