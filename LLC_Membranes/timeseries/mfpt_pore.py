#!/usr/bin/env python

""" Find the mean first passage time (MFPT) of a type of particle
"""

from LLC_Membranes.llclib import timeseries, file_rw, rand, stats
from LLC_Membranes.timeseries.fractional_levy_motion import FLM
from LLC_Membranes.timeseries import flm_sim_params
from LLC_Membranes.analysis.sfbm_parameters import SFBMParameters
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import tqdm
from multiprocessing import Pool
from scipy.sparse import csr_matrix as sparse_matrix
import warnings
from scipy.stats import norm
import levy
import fbm
import yaml
import argparse
import sys
import random


with warnings.catch_warnings():
    warnings.filterwarnings('ignore', "Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient.")
    warnings.filterwarnings('ignore', "SparseEfficiencyWarning")


def initialize():

    parser = argparse.ArgumentParser(description='Calculate mean first passage time from realizations of the specified'
                                                 'type of random walk.')

    parser.add_argument('-y', '--yaml', type=str, help='A .yaml configuration file. This is preferred over using'
                                                       'argparse. It leads to better reproducibility.')

    return parser


# TODO: should move these to llclib
class Hops:

    def __init__(self, distribution, params, n):
        """ Generate a sequence of hops, with a correlation structure of hurst parameter does not equal 0.5

        :param distribution: name of distribution. It should be defined in self.distributions
        :param params: parameters of the distribution in order (sigma, alpha, hurst, max_hop)
        :param n: number of points to generate when a call to hops is made (you will get best performance if a power of\
        2)

        :type distribution: str
        :type params: tuple or list
        :type n: int
        """

        self.name = distribution
        self.params = params
        self._distributions = {'gaussian': self._gaussian_hops, 'levy': self._levy_hops}
        self._set_params = {'gaussian': self._set_gaussian_params, 'levy': self._set_levy_params}
        self.length = n
        self.hops = self._distributions[self.name]
        self.set_params = self._set_params[self.name]

        if self.name == 'levy':

            print('Generating corrections to truncated Levy stable distribution...', end='', flush=True)
            self.hurst_correction = flm_sim_params.HurstCorrection()
            self.truncation_correction = flm_sim_params.TruncateLevy()
            corrected_max_hop = self.truncation_correction.interpolate(self.params[2], self.params[1] - .01,
                                                                       self.params[3], self.params[0])
            self.flm = FLM(self.hurst_correction.interpolate(self.params[2], self.params[1]), self.params[1],
                           scale=self.params[0],
                           N=self.length, M=4, correct_hurst=False, truncate=corrected_max_hop, correct_truncation=False)
            print('Done!')

    def _gaussian_hops(self):
        """ Draw hop lengths from a Gaussian distribution with a correlation structure defined by the hurst parameter.
        If the hurst parameter equals 0.5, it will be the same as Brownian motion

        :return: sequence of cumulative hop lengths
        """

        return self.params[0] * fbm.FBM(self.length, self.params[2], method="daviesharte").fbm() \
               / ((1.0 / self.length) ** self.params[2])

    def _set_gaussian_params(self, params):

        params[2] = 0.5  # for speed while testing

        self.params = (params[1], self.params[1], params[2], self.params[3])

    def _levy_hops(self):
        """ Draw a sequence of hop lengths from a Levy stable distribution

        :return:
        """

        self.flm.generate_realizations(1, progress=False)

        return self.flm.realizations[0, :]

    def _set_levy_params(self, params):

        pass


class Dwell:

    def __init__(self, distribution, params):
        """ Generate a sequence of dwell times by drawing from the appropriate dwell time distibution

        :param distribution: name of dwell distribution. If None, there are no dwells and a hop will occur every step
        :param params: parameters of chosen dwell time distribution. Passed in order (alpha, lambda, lower_limit)
        """

        self.name = distribution
        self._distributions = {None: self._no_trapping, 'power': self._powerlaw,
                               'power_cut': self._powerlaw_cut, 'exponential': self._exponential}
        self._set_params = {None: self._set_no_trapping_params, 'power': self._set_powerlaw_params,
                            'power_cut': self._set_powerlaw_cut_params, 'exponential': self._set_exponential_params}

        self.dwells = self._distributions[self.name]

        self.params = params

        self.set_params = self._set_params[self.name]

    @staticmethod
    def _no_trapping(n):
        """ Return a bunch of ones since there is no trapping occuring

        :param n: number of samples to draw
        :type n: int
        """

        return np.ones(n)

    @staticmethod
    def _set_no_trapping_params(params):

        pass

    def _powerlaw(self, n):
        """ Draw dwell times from pure power law. Note that this has an infinite variance

        :param n: number of samples to draw
        :type n: int
        """

        return rand.random_powerlaw(self.params[0], ll=self.params[2], size=n, limit=None,
                                    discrete=False, exact=False)

    def _set_powerlaw_params(self, params):

        self.params = (1 + params[0], self.params[1], self.params[2])

    def _powerlaw_cut(self, n):
        """ Draw dwell times from a power law with an exponential cut-off

        :param n: number of samples to draw
        :type n: int
        """

        return rand.random_powerlaw_cutoff(self.params[0], self.params[1], xmin=self.params[2], size=n)

    def _set_powerlaw_cut_params(self, params):

        self.params = (1 + params[0], params[1], self.params[2])

    def _exponential(self, n):
        """ Draw dwell times from an exponential distribution

        :param n: number of samples to draw
        :type n: int
        """

        return rand.random_exponential(self.params[1], size=n, xmin=self.params[1])

    def _set_exponential_params(self, params):

        self.params = (self.params[0], params[0], self.params[2])


class MeanFirstPassageTime:

    def __init__(self, L, n, nT, pickled_parameters=None, hop_params=None, dwell_params=None, dt=1., nbins=25,
                 nt=1, save=True, concentration_profile=False):
        """ Generate trajectories from which to calculate mean first passage times. One can also use them to calculate
        the concentration profile within the pore

        Only 1 mode anomalous diffusion trajectories are implemented

        :param L: length of 1D path (a pore for example) that particle follows (nm).
        :param n: number of particle trajectories to generate. It should be at least the number expected to be in a \
        pore at steady state, but more is always better.
        :param nT: number of steps more particle trajectory. Make this a multiple of 2 for the best performance. Note \
        that there are better ways to do this for pure Brownian motion, but this script is built with correlation in \
        mind.
        :param pickled_parameters: load hop and dwell paramters from pickle file generated by SFBMParameters class. This \
        will override hop_params and dwell_params
        :param hop_params: dictionary with keys distribution, sigma, alpha and hurst
        :param dwell_params: dictionary with keys distribution and alpha
        :param dt: timestep to use while contsructing trajectories
        :param sigma: width of hop distribution per unit time. The width of the hop distribution for the draws made each
        timestep is sqrt(dt) * sigma
        :param nbins: number of bins in concentration profile
        :param nt: number of threads to use when generating trajectories
        :param save: save trajectories to disk
        :param load: load trajectories saved to disk
        :param concentration_profile: set to True if you want to plot the concentration profile at the end. This
        requires all trajectories to be saved.

        :type L: float
        :type n: int
        :type hop_params: dict
        :type dwell_params: dict
        :type dt: float
        :type sigma: float
        :type nbins: int
        :type nt: int
        :type save: bool
        :type load: bool
        :type concentration_profile: bool
        """

        # trajectory building parameters
        self.length = L
        self.ntraj = n
        self.trajectories = []
        self.passage_times = []
        self.time = []
        self.dt = dt
        self.passage_times = []
        self.nt = nt
        self.discretization_points = discretization_points
        self.nbins = nbins
        self.resample = False

        self.store_traj = False
        if concentration_profile:
            self.store_traj = True

        hop_dist = hop_params['distribution']
        dwell_dist = dwell_params['distribution']

        # conform names of distribution
        hop_dist_name = None
        if hop_dist.lower() in ['gaussian', 'normal', 'brownian']:
            hop_dist_name = 'gaussian'
        elif hop_dist.lower() in ['levy']:
            hop_dist_name = 'levy'

        dwell_dist_name = None
        if dwell_dist is None:
            dwell_dist_name = None
        elif dwell_dist.lower() in ['power', 'powerlaw']:
            dwell_dist_name = 'power'
        elif dwell_dist.lower() in ['power_cut', 'powercut']:
            dwell_dist_name = 'power_cut'
        elif dwell_dist.lower() in ['exponential', 'exp']:
            dwell_dist_name = 'exponential'

        hop_p = None
        dwell_p = None

        if pickled_parameters is not None:
            # This can be improved but first requires improvements to SFBMParameters

            print('Loading Parameters from %s...' % pickled_parameters, end='', flush=True)
            measured_params = file_rw.load_object(pickled_parameters)
            print('Done!')

            self.hop_parameter_distribution = measured_params.hop_parameters[0]  # [0] is index of mode (only 1 mode works)
            self.dwell_parameter_distribution = measured_params.dwell_parameters[0]
            self.hurst_distribution = measured_params.hurst_distribution[0]

            if hop_dist_name == 'gaussian':

                sig = random.choice(self.hop_parameter_distribution)[1]
                hurst = random.choice(self.hurst_distribution)
                #hurst = 0.5

                hop_p = (sig, 2, hurst, measured_params.max_hop)

            elif hop_dist_name == 'levy':

                p = random.choice(self.hop_parameter_distribution)

                hop_p = (p[2], 2, p[0], measured_params.max_hop)

            if dwell_dist_name == 'power':

                dwell_alpha = random.choice(self.dwell_parameter_distribution)

                dwell_p = (1 + dwell_alpha, 0, measured_params.dwell_lower_limit[0])

            elif dwell_dist_name == 'power_cut':

                p = random.choice(self.dwell_parameter_distribution)

                dwell_p = (1 + p[0], p[1], measured_params.dwell_lower_limit[0])

            self.resample = True

        else:

            hop_p = (np.sqrt(self.dt) * hop_params['sigma'], hop_params['alpha'], hop_params['hurst'], hop_params['max_hop'])
            dwell_p = (dwell_params['alpha'], dwell_params['lambda'], dwell_params['lower_limit'])

        # print(dwell_dist_name)
        # print(hop_p)
        # print(dwell_p)
        # exit()

        self.hop_dist = Hops(hop_dist_name, hop_p, nT)
        self.dwell_dist = Dwell(dwell_dist_name, dwell_p)

        print('Generating Trajectories...')
        self.generate_trajectories(save=save)

        # flux simulation arrays
        self.positions = None
        self.time_uniform = None
        self.flux_in = None
        self.flux_out = None
        self.pore_concentration = 0  # concentration at pore entrance
        self.dz = 0
        self.steps = 0
        self.concentration = np.zeros([0, self.nbins])
        self.nparticles_in_pore = np.array([])

        # for animated plotting
        self.patches = None
        self.bins = None
        self.ax = None
        self.fig = None

    def _trajectory_realizations(self, ntraj):
        """ Generate realizations of random walk with predefined hop length and dwell time distributions

        :param ntraj: number of realizations to generate.
        :type ntraj: int
        """

        np.random.seed()  # need to set a different random seed for each thread

        trajectories = []
        times = []
        passage_times = []

        n = 0  # a diagnostic variable to track how many trajectories are generated in total

        while len(passage_times) < ntraj:

            if self.resample:

                hops = random.choice(self.hop_parameter_distribution)
                dwells = random.choice(self.dwell_parameter_distribution)
                hurst = random.choice(self.hurst_distribution)

                self.hop_dist.set_params(hops + [hurst])
                self.dwell_dist.set_params(dwells)

            walk = self.hop_dist.hops()

            # find when the cross-over occurs
            cross = timeseries.switch_points(np.abs(walk) >= self.length)[1] + 2

            walk = walk[:cross]

            traj = self._crossings(walk)

            if self.store_traj:

                for i, t in enumerate(traj):
                    trajectories.append(t)
                    time = self.dt * np.cumsum(self.dwell_dist.dwells(len(t)))
                    times.append(time - time[0])

                if np.abs(walk).max() >= self.length:
                    passage_times.append(times[-1][-1])

            else:

                if np.abs(walk).max() >= self.length:
                    time = np.cumsum(self.dwell_dist.dwells(len(traj[-1])))
                    passage_times.append(time[-1] * self.dt)

            n += 1

            print('\r%d/%d trajectories' % (len(passage_times), ntraj), end='')

        return trajectories, times, passage_times, n

    @staticmethod
    def _crossings(x):
        """ Identify whenever a

        :param x:
        :return:
        """

        start = 0
        end = 0
        traj = []

        while start < (x.size - 1):

            if x[start] * x[start + 1] < 0:
                end = start + 1
            else:
                end += (timeseries.switch_points(x[start:] > 0)[1] + 2)

            traj.append(np.abs(x[start:end]))
            start = end

        return traj

    def generate_trajectories(self, save=True, savename='trajectories.npz'):
        """ Generate Brownian trajectories

        :param nt: number of threads
        """

        pool = Pool(processes=self.nt)

        trajectories_per_thread = self.ntraj // self.nt

        arguments = [trajectories_per_thread for _ in range(self.nt)]

        result = pool.map(self._trajectory_realizations, arguments)

        pool.close()
        pool.join()

        n = 0
        for thread in range(self.nt):
            self.trajectories += result[thread][0]
            self.time += result[thread][1]
            self.passage_times += result[thread][2]
            n += result[thread][3]

        percent_crossed_pore = 100 * trajectories_per_thread * self.nt / n

        self.passage_times = np.array(self.passage_times)

        print('\n%.2f %% of trajectories reached L' % percent_crossed_pore)

        if save:
            file_rw.save_object((self.trajectories, self.time), 'trajectories.pl')

    def mfpt(self, nboot, percentile=99.5):
        """ Bootstrap the distribution of mean first passage times and report where the maximum occurs

        :param nboot: number of bootstrap trials
        :param percentile: plot all passage times up to this percentile of the distribution
        
        :type nboot: int
        :type percentile: float
        """

        max_pt = np.zeros(nboot)
        limits = (self.passage_times.min(), np.percentile(self.passage_times, percentile))
        bar_heights = np.zeros(self.nbins)

        for b in range(nboot):

            ptimes = np.random.choice(self.passage_times, size=self.passage_times.size, replace=True)

            n, edges = np.histogram(ptimes, bins=self.nbins, range=limits, density=True)
            bin_width = edges[1] - edges[0]
            bin_centers = [x + bin_width / 2 for x in edges[:-1]]

            max_pt[b] = bin_centers[np.argmax(n)]
            bar_heights += n

        bar_heights /= nboot

        plt.bar(bin_centers, bar_heights, bin_width)
        # from scipy.signal import savgol_filter as sf
        # plt.plot(bin_centers, sf(bar_heights, 19, 5), '--', lw=2, color='black')

        print('Mean first passage time: %.2f +/- %.2f ns' % (max_pt.mean(), max_pt.std()))

        plt.ylabel('Probability', fontsize=14)
        plt.xlabel('Passage time (ns)', fontsize=14)
        plt.tick_params(labelsize=14)
        plt.tight_layout()
        # self.mfpt_ecdf()

        plt.savefig('mfpt_L%s.pdf' % self.length)
        #plt.show()

    def mfpt_ecdf(self):
        """ Calculate the emperical cumulative distribution function of the passage times.
        """

        cdf = stats.Cdf(self.passage_times)

        # plt.figure()
        # plt.plot(cdf.xs, cdf.ys)

        x, uniq_ndx = np.unique(cdf.xs, return_index=True)
        y = cdf.ys[uniq_ndx]

        from scipy.interpolate import PchipInterpolator as pchip
        from scipy.signal import savgol_filter as sf

        smooth = sf(y, 201, 3)
        # plt.plot(x, smooth, '--', lw=2, color='black')
        #
        # #spline = pchip(x[::100], y[::100])
        #
        # plt.plot(x, y)
        # plt.show()
        # exit()
        # xinterp = np.linspace(x[0], x[-1], 500)
        # yinterp = spline(xinterp)

        # plt.plot(x, cdf.ys[uniq_ndx])
        # plt.show()
        # exit()

        #plt.figure()
        dx = x[1:] - x[:-1]
        deriv = np.diff(smooth) / dx
        file_rw.save_object(self.passage_times, 'ptimes.pl')
        plt.plot(x, sf(y, 501, 3, deriv=1), '--', label='Approximate Analytical PDF', lw=2, color='black')
        plt.show()
        exit()

    def save_passage_times(self, name='passage_times.pl'):

        file_rw.save_object(self.passage_times, name)

    def discretize(self, buffer=1):
        """ Discretize timeseries into equal length segments. Dwell time distributions are continuous. This is only
        useful for flux calculations because the concentration profile doesn't depend on the dwell time distribution.

        :param buffer: to discretize, the time series will be divided into equal length segments. This parameter is \
        multiplied by the number of observed time points in order to determine the number of segments. A higher value \
        will give a higher resolution trajectories with respect to non-uniform dwell times.

        :type buffer: int
        """

        longest = max([len(x) for x in self.trajectories])
        longest_time = max([x[-1] for x in self.time])
        self.time_uniform = np.linspace(0, longest_time, longest * buffer)
        self._discretize_time()

    def _discretize_time(self):
        # TODO: increase bins in time_uniform (e.g. make time steps one tenth of dt)

        ntraj = len(self.trajectories)
        nperthread = ntraj // self.nt

        arguments = []
        for n in range(self.nt):
            if n < (self.nt - 1):
                arguments.append((n * nperthread, (n + 1) * nperthread))
            else:
                arguments.append((n * nperthread, (n + 1) * nperthread + (ntraj % self.nt)))

        pool = Pool(processes=self.nt)
        pool.starmap(self._discretize_group, arguments)
        pool.close()

    def _discretize_group(self, start, end, max_memory=4 * 10 ** 9):

        for i in tqdm.tqdm(range(start, end)):

            #t = self.trajectories[i]
            t = self.time[i]

            last = np.argmin(np.abs(t[-1] - self.time_uniform))
            tu = self.time_uniform[:(last + 1)]

            # find the indices of the time point closest to the interpolated uniform time series
            if 8 * t.size * tu.size < max_memory:
                time_index = np.argmin(np.abs(t[:, np.newaxis] - tu), axis=0)
            else:
                # slower but less issues with memory
                time_index = []
                for x in tu:
                    time_ndx = np.argmin(np.abs(x - t))
                    time_index.append(time_ndx)

                time_index = np.array(time_index)

            # make sure that jumps don't occur in interpolated trajectory until jumps occur in real trajectory
            time_index[np.where(tu - t[time_index] < 0)[0]] -= 1
            time_index[np.where(time_index < 0)] += 1  # Relax above restriction for first time point

            # self.visualize_discretization(self.trajectories[i], t)
            # self.visualize_discretization(self.trajectories[i][time_index], tu)
            # plt.show()
            # exit()

            self.trajectories[i] = self.trajectories[i][time_index]
            self.time[i] = tu

    @staticmethod
    def visualize_discretization(z, time):

        # for visualizing hops
        trajectory_hops = np.zeros([2 * len(time) - 1, 2])

        trajectory_hops[1::2, 0] = time[1:]
        trajectory_hops[2::2, 0] = time[1:]

        trajectory_hops[::2, 1] = z
        trajectory_hops[1:-1:2, 1] = z[:-1]
        trajectory_hops[-1, 1] = z[-1]
        plt.plot(trajectory_hops[:, 0], trajectory_hops[:, 1])

    def concentration_from_histogram(self):
        """ Generate a concentration profile by histogramming the atomic """

        data = []
        for t in self.trajectories:
            data += t.tolist()

        #file_rw.save_object(data, 'data.pl')

        n, edges = np.histogram(data, bins=self.nbins, range=(0, self.length), density=True)
        bin_width = edges[1] - edges[0]
        bin_centers = [x + bin_width / 2 for x in edges[:-1]]

        plt.figure()
        plt.plot(bin_centers, n, lw=2)
        plt.fill_between(bin_centers, np.zeros(n.size), n, alpha=0.5)
        #plt.hist(data, self.nbins, range=(0, self.length), density=True)
        plt.xlabel('Pore coordinate (nm)', fontsize=14)
        plt.ylabel('Probability', fontsize=14)
        plt.tick_params(labelsize=14)
        plt.tight_layout()
        plt.show()

    def simulate_flux(self, pore_concentration=1, dz=0.1, steps=2000, measure_flux=False):
        """ Simulate a flux into the pore, starting new particles after a specified time
        """

        self.pore_concentration = pore_concentration

        longest = max([len(x) for x in self.trajectories])
        total_trajectories = len(self.trajectories)
        print('%d total trajectories generated' % total_trajectories)

        # discretize position versus time
        self.positions = sparse_matrix((self.ntraj, longest + steps * self.pore_concentration))

        self.dz = dz

        if self.dwell_dist.name is not None:
            print('Discretizing Trajectory...')
            n = 1
            self.time_uniform = np.arange(longest * n) * self.dt / n
            self.discretize_time()

        self.concentration_from_histogram()
        exit()

        self.steps = steps

        if measure_flux:
            self.flux_out = np.zeros(self.steps)

        traj_no = 0
        start = 0
        naddtot = []
        inlet = np.zeros(steps)
        dwell_times = []

        print('Simulating flux by holding the interface concentration constant at %d particles' % self.pore_concentration)
        for step in tqdm.tqdm(range(steps), unit=' Time Steps'):

            nadd = 0
            if step > 0:
                #print(inlet[:step], inlet[:step].mean())
                if inlet[:step].mean() < pore_concentration:
                    nadd = 1  # This will need to be changed to handle higher inlet concentration
            else:
                nadd = 1

            # figure out how many solutes to add in order to keep concentration constant
            #if self._get_inlet_concentration(start) == 0 and self._get_inlet_concentration(start + 1) == 0:
            #    nadd = 1
            #else:
            #    nadd = 0

            
            #c = self._get_inlet_concentration(start)
            #nadd = self.pore_concentration - c
            #inlet.append(c)

            #if nadd < 0:  # pore concentration too high. Don't add any
            #    nadd = 0

            naddtot.append(nadd)

            if traj_no + nadd > self.ntraj:
                print("Cleaning position matrix ...")

                previous_step = self.nparticles_in_pore.size
                # Record the total number of particles in the pore before modifying position matrix
                self._update_nparticles_in_pore(step - previous_step)

                # Record concentration profile in pore before modifying position matrix
                self._update_concentration_profile(step - previous_step)

                # Remove trajetories that are already finished
                self.positions, traj_no = self._clean_position_matrix(step - previous_step)
                start = 0

            for particle in range(traj_no, traj_no + nadd):

                random_trajectory_index = np.random.choice(total_trajectories)  # randomly choose a trajectory

                traj = self.trajectories[random_trajectory_index]
                if traj[-1] >= self.length:
                    dwell_times.append(len(traj))
                # print(traj)
                # exit()
                time = self.time[random_trajectory_index]
                last = np.argmin(np.abs(time[-1] - self.time_uniform))
                # tu = self.time_uniform[:(last + 1)]

                # find the indices of the time point closest to the interpolated uniform time series
                # time_index = np.argmin(np.abs(time[:, np.newaxis] - tu), axis=0)
                # # make sure that jumps don't occur in interpolated trajectory until jumps occur in real trajectory
                # time_index[np.where(tu - time[time_index] < 0)[0]] -= 1
                # time_index[np.where(time_index < 0)] += 1  # Relax above restriction for first time point
                #
                # interpolated = traj[time_index]  # TODO: this has issues for trajectories starting at zero

                #self.positions[particle, start:interpolated.size + start] = interpolated
                self.positions[particle, start:(traj.size + start)] = traj

                # print(self._get_inlet_concentration(start))

            if measure_flux:

                self.flux_out[step] = self._count_flux_out(start)

            inlet[step] = self._get_inlet_concentration(start)
            #print(inlet[step], nadd)
            #data = self.positions[:, start].data
            #print(data)
            #print(np.where(np.logical_and(data >=0, data < self.dz))[0].size)
            #if step == 10:
            #    exit()

            traj_no += nadd
            start += 1

        print(naddtot[equil:])
        print(np.mean(naddtot[equil:]))

        print(sum(naddtot) / (steps * self.dt))
        plt.plot(naddtot)
        #plt.plot(self.flux_out)
        plt.figure()
        plt.plot(inlet)
        plt.title('Inlet Concentration')
        file_rw.save_object(inlet, 'inlet.pl')
        file_rw.save_object(dwell_times, 'dwell_times.pl')

        plt.figure()
        plt.hist(dwell_times, bins=25)

        previous_step = self.nparticles_in_pore.size
        self._update_nparticles_in_pore(steps - previous_step)
        self._update_concentration_profile(steps - previous_step)

        #self.animate_concentration_profile()

    def simulate_flux_constant_particles(self, n=50, steps=2000, measure_flux=False):
        """ Simulate a flux into the pore, starting new particles after a specified time
        """

        longest = max([len(x) for x in self.trajectories])
        total_trajectories = len(self.trajectories)

        # discretize position versus time
        self.positions = sparse_matrix((self.ntraj, longest + steps))

        self.time_uniform = np.arange(self.positions.shape[1]) * self.dt
        self.steps = steps

        if measure_flux:
            self.flux_out = np.zeros(self.steps)

        traj_no = 0
        start = 0
        naddtot = []

        dwell_times = []
        print('Simulating flux by holding the nubmer of particles constant at %d' % n)
        for step in tqdm.tqdm(range(steps), unit=' Time Steps'):

            if step > 0:
                nadd = n - self._get_nparticles_in_pore(start)
                if nadd < 0:
                    print(nadd)
                    exit()
            else:
                nadd = n

            naddtot.append(nadd)

            if traj_no + nadd > self.ntraj:
                print("Cleaning position matrix ...")

                previous_step = self.nparticles_in_pore.size
                # Record the total number of particles in the pore before modifying position matrix
                self._update_nparticles_in_pore(step - previous_step)

                # Record concentration profile in pore before modifying position matrix
                self._update_concentration_profile(step - previous_step)

                # Remove trajetories that are already finished
                self.positions, traj_no = self._clean_position_matrix(step - previous_step)
                start = 0

            for particle in range(traj_no, traj_no + nadd):

                random_trajectory_index = np.random.choice(total_trajectories)  # randomly choose a trajectory

                traj = self.trajectories[random_trajectory_index]
                if traj[-1] >= self.length:
                    dwell_times.append(len(traj))

                time = self.time[random_trajectory_index]
                last = np.argmin(np.abs(time[-1] - self.time_uniform))
                # tu = self.time_uniform[:(last + 1)]

                # find the indices of the time point closest to the interpolated uniform time series
                # time_index = np.argmin(np.abs(time[:, np.newaxis] - tu), axis=0)
                # # make sure that jumps don't occur in interpolated trajectory until jumps occur in real trajectory
                # time_index[np.where(tu - time[time_index] < 0)[0]] -= 1
                # time_index[np.where(time_index < 0)] += 1  # Relax above restriction for first time point
                #
                # interpolated = traj[time_index]  # TODO: this has issues for trajectories starting at zero

                # self.positions[particle, start:interpolated.size + start] = interpolated
                self.positions[particle, start:(traj.size + start)] = traj

                # print(self._get_inlet_concentration(start))

            if measure_flux:
                self.flux_out[step] = self._count_flux_out(start)

            traj_no += nadd
            start += 1

        file_rw.save_object(self.flux_out, 'flux.pl')

        plt.plot(naddtot)

        file_rw.save_object(dwell_times, 'dwell_times.pl')

        plt.figure()
        plt.hist(dwell_times, bins=25)

        previous_step = self.nparticles_in_pore.size
        self._update_nparticles_in_pore(steps - previous_step)
        self._update_concentration_profile(steps - previous_step)

    def _get_nparticles_in_pore(self, step):

        data = self.positions[:, step].data
        return np.where(np.logical_and(data >= 0, data < self.length))[0].size

    def _get_inlet_concentration(self, step):
        """ Number of particles in inlet slab

        :param step: step number

        :type step: int

        :return: number of particles in inlet slab
        """

        # this gets progressively slower as sparse matrix grows
        # return np.where(self.positions.tocsr()[:, step].data < self.dz)[0].size

        # for lil_matrix
        # return np.where(np.array(self.positions[:, step].T.data[0]) < self.dz)[0].size

        # when self.positions is already a csr matrix
        #return np.where(self.positions[:, step].data < self.dz)[0].size
        data = self.positions[:, step].data
        return np.where(np.logical_and(data >=0, data < self.dz))[0].size

    def _clean_position_matrix(self, step):
        """ Create a new position matrix, discarding trajectories that have left the pore already.
        """

        positions = sparse_matrix(self.positions.shape)
        nonzero = self.positions[:, step].nonzero()[0]

        positions[:nonzero.size, :(positions.shape[1] - step)] = self.positions[nonzero, step:]

        first_free_slot = nonzero.size

        return positions, first_free_slot

    def _count_flux_out(self, step):
        """ Return the number of particles which left the pore at x = L during this step
        """

        return np.where(self.positions[:, step].data >= self.length)[0].size
        #return np.nonzero(self.positions[:, step] >= self.length)[0].size

    def _update_nparticles_in_pore(self, step):
        """ Update the array keeping track of the number of particles in the pore

        :param step: the frame up until which data should analyzed

        :type step: int
        """

        nonzero = self.positions[:, :step].getnnz(axis=0)

        self.nparticles_in_pore = np.concatenate((self.nparticles_in_pore, nonzero))

    def _update_concentration_profile(self, step):

        concentration = np.zeros([step, self.nbins])
        # data = self.positions.T.data

        for t in range(step):
            concentration[t, :] = np.histogram(self.positions[:, t].data, self.nbins, range=(0, self.length), density=False)[0]

        self.concentration = np.concatenate((self.concentration, concentration))

    def plot_number_particles(self, show=False):

        plt.figure()

        plt.plot(self.time_uniform[:self.steps], self.nparticles_in_pore)

        if show:
            plt.show()

    def plot_average_concentration(self, equil=500, show=False, theoretical=True, save=True):
        """ Plot the average concentration profile along the pore

        :param equil: the frame at which the flux can be considered to be equilibrated

        :type equil: int
        """

        plt.figure()
        bar_edges = np.linspace(0, self.length, self.nbins + 1)

        bin_width = bar_edges[1] - bar_edges[0]
        bar_locations = [b + (bin_width / 2) for b in bar_edges[:-1]]
        plt.bar(bar_locations, self.concentration[equil:, :].mean(axis=0), bin_width, align='center')
        plt.xlabel('Distance along pore axis (nm)', fontsize=14)
        plt.ylabel('Concentration', fontsize=14)

        if theoretical:
            # theoretical profile based on measured flux  c = c0 * (1 - x / L)
            x = np.linspace(0, self.length, 1000)
            c = self.pore_concentration * (1 - x / self.length)
            plt.plot(x, c, '--', color='black', lw=2, label='Theoretical Profile')
            plt.legend(fontsize=14)

        plt.tick_params(labelsize=14)
        plt.tight_layout()
        #plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/stochastic_transport/supporting_figures/brownian_conc_profile.pdf')

        # x = np.linspace(0, self.length, 1000)
        # c0 = self.pore_concentration
        # j = self.pore_concentration
        # plt.plot(x, c0 - j * x)

        if show:
            plt.show()

    def animate_concentration_profile(self, bins=25, step=1):

        self.fig, self.ax = plt.subplots()
        nonzero = np.nonzero(self.positions[:, step] >= 0)[0]

        n, self.bins, self.patches = plt.hist(self.positions[nonzero, 0], bins, range=(0, self.length), density=True)

        ani = animation.FuncAnimation(self.fig, self._animate, blit=True,
                                      frames=iter(np.arange(0, self.positions.shape[1], step)), repeat=False)

        # ani.save('/home/bcoscia/brownian_impulse.gif', writer='imagemagick', fps=3)
        plt.show()

    def _animate(self, frame):

        try:
            nonzero = np.nonzero(self.positions[:, frame] >= 0)[0]
            n, _ = np.histogram(self.positions[nonzero, frame], self.bins, range=(0, self.length), density=True)
        except FloatingPointError:  # happens if there is no data
            n = np.zeros(self.bins)

        for rect, h in zip(self.patches, n):
            rect.set_height(h)

        self.ax.set_title('Frame %d' % frame)

        self.ax.relim()
        self.ax.autoscale_view()
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

        return self.patches


if __name__ == "__main__":

    """ NOTES
    - include first step of trajectory (at zero) for flux simulations. But exclude from histogram method
    - script now maintains concentration based on average of all previoius frames. Might want to update to a moving window since might only want to average equilibrium portion with equilibrium concentrations
    - pore concentration fixed to read 0 < x < dz
    - need to make dz independent from number of bins in concentration profile. Can be done by normalizing conc profile
    - larger dz makes it harder to maintain concentration (?)
    """

    args = initialize().parse_args()

    if args.yaml:
        with open(args.yaml, 'r') as yml:
            cfg = yaml.load(yml)
    else:
        sys.exit('You must provide a yaml configuration file.')

    L = cfg['pore_length']
    ntraj = cfg['ntraj']  # this is the number of trajectories that actually make it to the end
    sigma = cfg['hops']['sigma']
    dt = cfg['dt']
    nt = cfg['nthreads']
    hop_parameters = cfg['hops']
    dwell_parameters = cfg['dwells']
    nT = 2**cfg['nT']
    save_trajectories = cfg['save_trajectories']
    discretization_points = cfg['discretization_points']
    bins = cfg['bins']

    print('sigma per step: %.2f' % (sigma * np.sqrt(dt)))

    # flux stuff
    pore_conc = 1
    steps = 500000
    equil = int(steps/2)  # use 3/4 of the data
    dz = L / bins  # make dz independent of bins. Will just need to change normalization on plot
    print('dz: %.2f' % dz)
    nparticles = 5

    load = False

    mfpt = MeanFirstPassageTime(L, ntraj, nT, hop_params=hop_parameters, dwell_params=dwell_parameters, dt=dt,
                                nbins=bins, nt=nt, save=save_trajectories, pickled_parameters=cfg['pickled_parameters'],
                                concentration_profile=cfg['concentration'])
    mfpt.mfpt(cfg['nboot'])
    if cfg['concentration']:
        mfpt.concentration_from_histogram()

    if cfg['save_passage_times']:
        mfpt.save_passage_times(name=cfg['save_passage_times_name'])

    exit()

    mfpt.simulate_flux(pore_concentration=pore_conc, dz=dz, steps=steps, measure_flux=True)
    #mfpt.simulate_flux_constant_particles(n=nparticles, dt=dt, steps=steps, measure_flux=True)
    print(mfpt.flux_out[equil:].mean() / dt)

    mfpt.plot_average_concentration(equil=equil)
    mfpt.plot_number_particles(show=True)
    exit()

    plt.hist(mfpt.passage_times, bins=50, range=(0, 10000))

    # for t in range(3):
    #     plt.plot(mfpt.trajectory_hops[t, :, 0], mfpt.trajectory_hops[t, :, 1])
    plt.show()
