#!/usr/bin/env python

""" Find the mean first passage time (MFPT) of a type of particle
"""

from LLC_Membranes.llclib import timeseries, file_rw, rand, stats, fitting_functions
from LLC_Membranes.timeseries.fractional_levy_motion import FLM
from LLC_Membranes.timeseries import flm_sim_params
from LLC_Membranes.analysis.markov_state_dependent_dynamics import States, Chain
from LLC_Membranes.analysis.sfbm_parameters import SFBMParameters
from scipy.integrate import quad
from scipy.optimize import curve_fit
from scipy.stats import levy_stable
import matplotlib.pyplot as plt
import numpy as np
import tqdm
from multiprocessing import Pool
import warnings
import fbm
import yaml
import argparse
import sys
import random
import time as timer


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
            corrected_max_hop = self.truncation_correction.interpolate(self.params[0], self.params[1] - .01,
                                                                       self.params[2], self.params[3])
            self.flm = FLM(self.hurst_correction.interpolate(self.params[0], self.params[1]), self.params[1],
                           scale=self.params[3],
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

    def __init__(self, L, n, nT, model, pickled_parameters=None, hop_params=None, dwell_params=None, dt=1., nbins=25,
                 nt=1, save=True, concentration_profile=False, load_passage_times=None, timeout=None):
        """ Generate trajectories from which to calculate mean first passage times. One can also use them to calculate
        the concentration profile within the pore

        Only 1 mode anomalous diffusion trajectories are implemented

        :param L: length of 1D path (a pore for example) that particle follows (nm).
        :param n: number of particle trajectories to generate. It should be at least the number expected to be in a \
        pore at steady state, but more is always better.
        :param nT: number of steps more particle trajectory. Make this a multiple of 2 for the best performance. Note \
        that there are better ways to do this for pure Brownian motion, but this script is built with correlation in \
        mind.
        :param model: name of model to use in order to generate realizations (AD or MSDDM)
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
        :param concentration_profile: set to True if you want to plot the concentration profile at the end. This
        requires all trajectories to be saved.
        :param load_passage_times: load passage times saved to disk.

        :type L: float
        :type n: int
        :type nT: int
        :type model: str
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

        print('Simulating pore of length %.1f nm' % L)

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
        self.time_uniform = None
        self.model = model
        self.chains = None  # will hold Chain object if model is MSDDM
        self.hop_realization = None
        self.dwell_realization = None
        self.nT = nT
        self.bound = None
        self.timeout = timeout

        self.store_traj = False
        if concentration_profile:
            self.store_traj = True

        if load_passage_times is None:

            if self.model.lower() in ['ad', 'anomalous diffusion']:

                self._load_anomalous_diffusion_params(hop_params, dwell_params, pickled_parameters)
                self.hop_realization = self._ad_hops
                self.dwell_realization = self._ad_dwells

                print('Generating Trajectories...')

            elif self.model.lower() in ['msddm']:

                self._initialize_msddm(hop_params, pickled_parameters)
                self.hop_realization = self._msddm_hops
                self.dwell_realization = self._msddm_dwells

            else:

                raise Exception("The model '%s' is not implemented" % self.model)

            self.generate_trajectories(save=save)

        else:

            self.passage_times = file_rw.load_object(load_passage_times)

    def _initialize_msddm(self, hop_params, pickled_parameters, speedup=True):

        if pickled_parameters is None:

            raise Exception('You must provide the parameters of the MSDDM in a pickle file formatted as output'
                            'by markov_state_dependent_dynamics.py')

        states = file_rw.load_object(pickled_parameters)

        self.bound = hop_params['max_hop']

        fbm = False
        try:
            fbm = hop_params['fbm']
        except KeyError:
            pass

        self.chains = Chain(states.count_matrix, states.fit_params, hurst_parameters=states.hurst,
                            emission_function=levy_stable, speedup=speedup, bound=self.bound, fbm=fbm)

    def _load_anomalous_diffusion_params(self, hop_params, dwell_params, pickled_parameters):

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

            self.hop_parameter_distribution = measured_params.hop_parameters[
                0]  # [0] is index of mode (only 1 mode works)
            self.dwell_parameter_distribution = measured_params.dwell_parameters[0]
            self.hurst_distribution = measured_params.hurst_distribution[0]

            if hop_dist_name == 'gaussian':

                sig = random.choice(self.hop_parameter_distribution)[1]
                hurst = random.choice(self.hurst_distribution)
                # hurst = 0.5

                hop_p = (sig, 2, hurst, measured_params.max_hop)

            elif hop_dist_name == 'levy':

                p = random.choice(self.hop_parameter_distribution)
                hurst = random.choice(self.hurst_distribution)

                # H, alpha, max, scale
                hop_p = (hurst, p[0], measured_params.max_hop, p[2])

            if dwell_dist_name == 'power':

                dwell_alpha = random.choice(self.dwell_parameter_distribution)

                dwell_p = (1 + dwell_alpha, 0, measured_params.dwell_lower_limit[0])

            elif dwell_dist_name == 'power_cut':

                p = random.choice(self.dwell_parameter_distribution)

                dwell_p = (1 + p[0], p[1], measured_params.dwell_lower_limit[0])

            self.resample = True

        else:

            hop_p = (
            np.sqrt(self.dt) * hop_params['sigma'], hop_params['alpha'], hop_params['hurst'], hop_params['max_hop'])
            dwell_p = (dwell_params['alpha'], dwell_params['lambda'], dwell_params['lower_limit'])

        self.hop_dist = Hops(hop_dist_name, hop_p, self.nT)
        self.dwell_dist = Dwell(dwell_dist_name, dwell_p)

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
            #start = timer.time()
            walk = self.hop_realization()

            # find when the cross-over occurs
            cross = timeseries.switch_points(np.abs(walk) >= self.length)[1] + 2

            walk = walk[:cross]
            print(walk[-1])

            traj = self._crossings(walk)

            if self.store_traj:

                for i, t in enumerate(traj):

                    trajectories.append(t)

                    dwells = self.dwell_realization(len(t))

                    time = self.dt * np.cumsum(dwells)
                    times.append(time - time[0])

                if np.abs(walk).max() >= self.length:
                    passage_times.append(times[-1][-1])

            else:

                if np.abs(walk).max() >= self.length:

                    dwells = self.dwell_realization(len(traj[-1]))
                    # dwells = np.ones(len(traj[-1]))

                    time = np.cumsum(dwells)
                    #
                    # plt.plot(time, traj[-1])
                    # plt.show()
                    # exit()

                    passage_times.append(time[-1] * self.dt)

            n += 1

            print('\r%d/%d trajectories' % (len(passage_times), ntraj), end='', flush=True)
            #print(timer.time() - start)

        return trajectories, times, passage_times, n

    def _ad_hops(self):

        if self.resample:
            hops = random.choice(self.hop_parameter_distribution)
            dwells = random.choice(self.dwell_parameter_distribution)
            hurst = random.choice(self.hurst_distribution)

            self.hop_dist.set_params(hops + [hurst])
            self.dwell_dist.set_params(dwells)

        return self.hop_dist.hops()

    def _ad_dwells(self, n):

        return self.dwell_dist.dwells(n)

    def _msddm_hops(self):

        # some hard sets for now
        self.chains.generate_realizations(1, self.nT, bound=self.bound, m=256, Mlowerbound=8, nt=1, quiet=True,
                                          timeout=self.timeout)

        return self.chains.chains[:, 0]

    @staticmethod
    def _msddm_dwells(n):

        return np.ones(n)

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

    def mfpt(self, nboot, method='peak', percentile=99.5, tstart=0, show=False):
        """ Bootstrap the distribution of mean first passage times and report where the maximum occurs

        :param nboot: number of bootstrap trials
        :param method: method of calculating MFPT. 'peak': define the peak of the passage time distribution as the \
        mean; 'fit_tails': Fit a stretch exponential function to the tail of the distribution in order to calculate a \
        mean of an incomplete passage time distribution; 'analytical': An analytical distribution of first passage times
        based on the Marathon runner problem.
        :param percentile: plot all passage times up to this percentile of the distribution
        :param tstart: Needed for 'fit_tails' method. time after which to start fitting
        :param show: show the fit to the histogram

        
        :type nboot: int
        :type percentile: float
        :type method: str
        :type tstart: int
        :type show: bool
        """

        methods = {'peak': self._mfpt_peak, 'fit_tails': self._mfpt_fit_tails, 'analytical': self._mfpt_analytical}

        methods[method](nboot, percentile, tstart=tstart)

    def _mfpt_analytical(self, nboot, percentile, **kwargs):

        self.passage_times /= 1000

        ptimes = np.random.choice(self.passage_times, size=self.passage_times.size, replace=True)

        hist, edges = np.histogram(ptimes, self.nbins, density=True)  # convert to microseconds
        bin_width = edges[1] - edges[0]
        bin_centers = np.array([i + bin_width / 2 for i in edges[:-1]])

        # very important to have a good guess. Might need to pass these unless I calculate MSD
        #p0 = [self.length, self.passage_times.mean() / self.length, 2.5]
        p0 = [self.length, 1, 2.5]
        # print(p0)
        # t = np.linspace(0.1, self.passage_times.max(), 1000)
        # plt.plot(t, fitting_functions.continuum_passage_time_distribution(t, p0[0], p0[1], p0[2]), '--',
        #          color='black', lw=2)
        # plt.show()
        # exit()

        epsilon = 0.00001  # tolerance in parameter
        bounds = [(self.length - epsilon, 0, 0), (self.length + epsilon, np.inf, np.inf)]

        popt = curve_fit(fitting_functions.continuum_passage_time_distribution, bin_centers, hist, p0=p0,
                         bounds=bounds)[0]

        #print(popt)

        t = np.linspace(0.1, self.passage_times.max(), 1000)
        mfpt = quad(fitting_functions.continuum_ptime_distribution_expected_value, 0.1, np.inf,
                    args=(popt[0], popt[1], popt[2]))[0]
        #print(mfpt)

        # if show:
        #
        #     plt.plot(t, fitting_functions.continuum_passage_time_distribution(t, popt[0], popt[1], popt[2]), '--',
        #              color='black', lw=2)
        #
        #     plt.hist(ptimes, self.nbins, density=True)
        #     plt.show()

    def _mfpt_fit_tails(self, nboot, percentile, **kwargs):
        """ Fit a stretch exponential function to the tail of an incomplete passage time distribution

        :param nboot: number of bootstrap trials
        :param percentile: plot up to this percentile of distribution.

        :type nboot: int
        :type percentile: float
        """

        # fit to data from tstart to the passage time distribution maximum
        tstart = kwargs['tstart']

        boot = np.zeros([nboot])
        for b in range(nboot):

            ptimes = np.random.choice(self.passage_times, size=self.passage_times.size, replace=True)

            # histogram data
            hist, edges = np.histogram(ptimes, self.nbins, density=False)  # convert to microseconds
            bin_width = edges[1] - edges[0]
            bin_centers = np.array([i + bin_width / 2 for i in edges[:-1]])
            ndx = np.where(bin_centers >= tstart)[0]  # bins from t start on

            p0 = [hist[ndx[0]], .5]  # initial guess at stretched exponential parameters
            bounds = [(0, 0), (np.inf, 1)]

            popt, _ = curve_fit(fitting_functions.stretched_exponential, bin_centers[ndx], hist[ndx], p0=p0,
                                bounds=bounds)

            area_front = sum([bin_width*hist[i] for i in range(ndx[0])])
            mean_front = ptimes[ptimes < bin_centers[-1]].mean()

            # numerically integrate. There is an exact form as well
            area_tail = quad(fitting_functions.stretched_exponential, bin_centers[-1], np.inf,
                             args=(popt[0], popt[1]))[0]
            mean_tail = quad(fitting_functions.stretched_exponential_expected_value, bin_centers[-1], np.inf,
                             args=(popt[0], popt[1], area_tail))[0]

            total_area = area_front + area_tail

            boot[b] = mean_front * (area_front / total_area) + mean_tail * (area_tail / total_area)

        print('Mean First Passage Time: %.2f +/- %.2f microseconds' % (boot.mean(), boot.std()))

        x = np.linspace(bin_centers[ndx[0]], bin_centers[ndx[-1]], 100)
        plt.plot(x, popt[0] * np.exp(-x**popt[1]), '--', color='black')
        plt.hist(self.passage_times, self.nbins, density=False)
        plt.show()

    def _mfpt_peak(self, nboot, percentile):

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
    model = cfg['model']
    timeout = cfg['timeout']

    print('sigma per step: %.2f' % (sigma * np.sqrt(dt)))

    # flux stuff
    pore_conc = 1
    steps = 500000
    equil = int(steps/2)  # use 3/4 of the data
    dz = L / bins  # make dz independent of bins. Will just need to change normalization on plot
    print('dz: %.2f' % dz)
    nparticles = 5

    load = None
    if cfg['load_passage_times']:
        load = cfg['save_passage_times_name']

    mfpt = MeanFirstPassageTime(L, ntraj, nT, model, hop_params=hop_parameters, dwell_params=dwell_parameters, dt=dt,
                                nbins=bins, nt=nt, save=save_trajectories, pickled_parameters=cfg['pickled_parameters'],
                                concentration_profile=cfg['concentration'], load_passage_times=load, timeout=timeout)

    #mfpt.mfpt(cfg['nboot'], method=cfg['mfpt']['method'], tstart=cfg['mfpt']['tstart'])

    if cfg['concentration']:
        mfpt.concentration_from_histogram()

    if cfg['save_passage_times']:
        mfpt.save_passage_times(name=cfg['save_passage_times_name'])
