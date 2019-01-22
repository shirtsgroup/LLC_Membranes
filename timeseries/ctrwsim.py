#!/usr/bin/env python

import argparse
import numpy as np
import tqdm
import sys
from LLC_Membranes.llclib import timeseries, fitting_functions, stats, sampling
from LLC_Membranes.analysis import Poly_fit
import matplotlib.pyplot as plt
from multiprocessing import Pool
import time as timer
from fbm_fork import fbm


def initialize():

    parser = argparse.ArgumentParser(description='Simulate fractional brownian motion and calculate its MSD')

    parser.add_argument('-hop', '--hop_length_distribution', default='gaussian', help='Functional form of hop length'
                                                                                      'distribution')
    # if more distributions are included, this will need to be more complicated depending what parameters are needed
    parser.add_argument('-hs', '--hop_sigma', default=1, type=float, help='Standard deviation of gaussian distribution '
                                                                          'used for drawing hop lengths')
    parser.add_argument('-H', '--hurst', default=0.5, type=float, help='Hurst parameter, used if hops are modeled with '
                                                                       'fractional brownian motion')
    parser.add_argument('-n', '--noise', default=0, type=float, help='Magnitude of gaussian noise to add to generated'
                                                                     'time series')
    parser.add_argument('-dwell', '--dwell_time_distribution', default='power', help='Functional form of dwell time'
                        'distribution (options: power, exponential')
    parser.add_argument('-steps', '--steps', default=1000, type=int, help='Number of steps to take for each'
                                                                          'independent trajectory')
    parser.add_argument('-ntraj', '--ntraj', default=100, type=int, help='Number of independent ctrw trajectories.')
    parser.add_argument('-ensemble', '--ensemble', action="store_true", help='Calculate MSD as ensemble average')
    parser.add_argument('-power_law', '--fit_power_law', action="store_true", help='Fit MSD to a power law')
    parser.add_argument('-linear', '--fit_line', action="store_true", help='Fit a line to the MSD')
    parser.add_argument('-alpha', '--alpha', default=0.5, type=float, help='Anomalous exponent')
    parser.add_argument('-lambda', '--lambda', default=0.5, type=float, help='Exponential decay rate')
    parser.add_argument('-dt', '--dt', default=1, type=float, help='Discrete time step for fixed length simulations')
    parser.add_argument('-fix_time', '--fix_time', action="store_true", help='Fix the total time of simulated '
                        'trajectories. Total length will be steps*dt')
    parser.add_argument('-b', '--nboot', default=200, help='Number of bootstrap trials')
    parser.add_argument('-ul', '--upper_limit', default=None, type=int, help='Upper limit on dwell-time length')
    parser.add_argument('-acf', '--autocorrelation', action="store_true", help='Plot step autocorrelation function')

    # parallelization
    parser.add_argument('-nt', '--nthreads', default=0, type=int, help='Number of threads to use for parallelized '
                                                                       'portions of the code.')

    return parser


class CTRW(object):

    def __init__(self, length, ntraj, hop_dist='gaussian', dwell_dist='power', hop_sigma=1, alpha=0.5, lamb=0.5,
                 padding=10, dt=1, nt=1, H=0.5):
        """ Initialize simulation of a continuous time random walk

        :param length: length of each simulated trajectory. If you fix the number of steps, this equals the number of
        steps. If you fix the time, the total length of the simulation is length * dt
        :param ntraj: number of independent trajectories to generate
        :param hop_dist: Method used to generate random hop lengths. "gaussian" draws randomly from a gaussian
        distribution, while "fbm" generates correlated hops as in fractional brownian motion.
        :param dwell_dist: Name of probability distribution function used to generate random dwell times
        :param hop_sigma: Sigma for Gaussian hop_dist random draws
        :param alpha: Anomalous exponent for power law random draws
        :param lamb: rate of decay for exponential random draws
        :param padding: multiplies number of discrete time points used to interpolate trajectories
        :param dt: time step for fixed time simulations
        :param nt: number of threads to use where parallelized
        :param H: Hurst parameter for fractional brownian motion. Equals 2*alpha for pure FBM. For H < 0.5, trajectories
        are subdiffusive, when H = 0.5, brownian motion is recovered.

        :type length: int
        :type ntraj: int
        :type hop_dist: str
        :type dwell_dist: str
        :type hop_sigma: float
        :type alpha: float
        :type lamb: float
        :type padding: int
        :type dt: float
        :type nt: int
        :type H: float
        """

        self.nsteps = length
        self.ntraj = ntraj
        self.hop_distribution = hop_dist
        self.hop_sigma = hop_sigma
        self.dwell_distribution = dwell_dist
        self.lamb = lamb
        self.alpha = alpha
        self.padding = padding
        self.dt = dt
        self.nt = nt
        self.H = H

        self.trajectories = np.zeros([self.ntraj, self.nsteps, 2])
        self.trajectory_hops = np.zeros([self.ntraj, 2 * self.nsteps - 1, 2])  # for visualization
        self.time_uniform = None
        self.z_interpolated = np.zeros([self.ntraj, self.nsteps*self.padding])  # separate from time_uniform to save memory

        self.msd = None
        self.fit_parameters = None
        self.bootstraps = None
        self.fit_cut = 1
        self.fit_start = 0
        self.acf = None  # autocorrelation function
        self.final_msd = None

        # Initialize multi-threading
        self.pbar = None

    def generate_trajectories(self, fixed_time=False, noise=0, ll=0.1, limit=None, distributions=None, discrete=False):
        """ Create trajectories by randomly drawing from dwell and hop distributions

        :param fixed_time: Propagate each trajectory until a certain wall time is reached
        :param noise: add gaussian noise to final trajectories
        :param ll: lower limit of power law distribution
        :param limit: upper limit on dwell time length (As implemented, this is not mathematically sound, so it's only
        for demonstration purposes)
        :param distributions: distributions of alpha and sigma values for dwell and hop length distributions
        respectively. Passed as 2-tuple of arrays where each array contains a possible values of each parameter.
        :param discrete: pull from discrete dwell time probability distributions

        :type fixed_time: bool
        :type noise: float
        :type ll: float
        :type limit: float
        :type distributions: tuple
        :type discrete: bool
        """

        if fixed_time:
            self.fixed_time_trajectories(ll=ll, distributions=distributions, discrete=discrete, noise=noise)
        else:
            self.fixed_steps_trajectories(noise=noise, nt=self.nt, ll=ll, limit=limit)

    def fixed_time_trajectories(self, ll=1, distributions=None, discrete=False, noise=0):
        """ Create trajectories of fixed length where dwell times and hop lengths are drawn from the appropriate
        distributions.

        :param ll: lower limit of power law distribution
        :param distributions: distributions of alpha and sigma values for dwell and hop length distributions
        respectively. Passed as 2-tuple of arrays where each array contains a possible values of each parameter.
        :param discrete: pull from discrete dwell time probability distributions

        :type ll: float
        :type distributions: tuple
        :type discrete: bool
        """

        self.time_uniform = np.linspace(0, self.nsteps, self.nsteps * self.padding)

        for t in tqdm.tqdm(range(self.ntraj)):

            if distributions is not None:

                if self.dwell_distribution == 'exponential':
                    self.lamb = np.random.choice(distributions[0])
                elif self.dwell_distribution == 'power':
                    self.alpha = np.random.choice(distributions[0])

                self.hop_sigma = np.random.choice(distributions[1])
                self.H = np.random.choice(distributions[2])

            time = [0]
            total_time = 0  # saves a lot of time

            while total_time < self.nsteps:

                # hop at random time intervals according to one of the following PDFs
                if self.dwell_distribution == 'exponential':
                    time.append(sampling.random_exponential_dwell(self.lamb))
                elif self.dwell_distribution == 'power':
                    time.append(sampling.random_power_law_dwell(1 + self.alpha, ll=ll, discrete=discrete))
                else:
                    sys.exit('Please enter a valid dwell time probability distribution')
                total_time += time[-1]

            time = np.cumsum(time)

            if self.hop_distribution in ['gaussian', 'Gaussian']:

                z = np.cumsum(np.random.normal(loc=0, scale=self.hop_sigma, size=len(time)))
                z -= z[0]  # untested

            elif self.hop_distribution in ['fbm', 'fractional', 'fraction_brownian_motion']:
                z = fbm.FBM(len(time), self.H, method="daviesharte", scale=self.hop_sigma).fbm()[:-1]  # automatically inserts zero at beginning of array

            else:
                sys.exit('Please enter a valid hop distance probability distribution')

            # for visualizing hops
            # trajectory_hops = np.zeros([2 * len(time) - 1, 2])
            #
            # trajectory_hops[1::2, 0] = time[1:]
            # trajectory_hops[2::2, 0] = time[1:]
            #
            # trajectory_hops[::2, 1] = z
            # trajectory_hops[1:-1:2, 1] = z[:-1]
            # trajectory_hops[-1, 1] = z[-1]
            # plt.plot(trajectory_hops[:, 0], trajectory_hops[:, 1])
            # plt.show()
            # exit()

            # make uniform time intervals with the same interval for each simulated trajectory
            self.z_interpolated[t, :] = z[np.digitize(self.time_uniform, time, right=False) - 1]

        self.time_uniform *= self.dt
        # plt.plot(trajectory_hops[:, 0]*self.dt, trajectory_hops[:, 1])
        # plt.plot(self.time_uniform, self.z_interpolated[-1, :])
        # plt.show()
        # exit()

    def fixed_steps_trajectories(self, noise=0, nt=1, ll=0.1, limit=None):
        """ Generate CTRW trajectories using a fixed number of steps

        :param noise: magnitude of gaussian noise to add to each trajectory
        :param nt: number of threads to use where parallelized
        :param ll: lower limit on dwell time length
        :param limit: upper limit on dwell time length (As implemented, this is not mathematically sound, so it's only
        for demonstration purposes)

        :type noise: float
        :type nt: int
        :type ll: float
        :type limit: float
        """

        print('Generating Trajectories...')
        for i in tqdm.tqdm(range(self.ntraj)):

            if self.hop_distribution == 'gaussian' or self.hop_distribution == 'Gaussian':
                z_position = np.cumsum(
                    np.random.normal(loc=0, scale=self.hop_sigma, size=self.nsteps))  # accumulate gaussian steps
            else:
                sys.exit('Please enter a valid hop distance probability distribution')

            self.trajectories[i, :, 1] = z_position - z_position[0]  # make initial z equal to 0

            # hop at random time intervals according to one of the following PDFs
            if self.dwell_distribution == 'exponential':
                time = sampling.random_exponential_dwell(self.lamb, size=self.nsteps)
            elif self.dwell_distribution == 'power':
                time = sampling.random_power_law_dwell(1 + self.alpha, size=self.nsteps, ll=ll, limit=limit)
            else:
                sys.exit('Please enter a valid dwell time probability distribution')

            time = np.cumsum(time)  # accumulate dwell times
            time -= time[0]

            self.trajectories[i, :, 0] = time

            # Add to array with all corners of hop distribution for visualization purposes
            self.trajectory_hops[i, 1::2, 0] = time[1:]
            self.trajectory_hops[i, 2::2, 0] = time[1:]

            self.trajectory_hops[i, ::2, 1] = self.trajectories[i, :, 1]
            self.trajectory_hops[i, 1:-1:2, 1] = self.trajectories[i, :-1, 1]
            self.trajectory_hops[i, -1, 1] = self.trajectories[i, -1, 1]

        print('Interpolating Trajectories...')
        # make uniform time intervals with the same interval for each simulated trajectory
        max_time = np.min(self.trajectories[:, -1, 0])
        self.time_uniform = np.linspace(0, max_time, self.nsteps*10)

        if nt > 1:
            # self.pbar = tqdm.tqdm(total=self.ntraj)
            pool = Pool(nt)
            for i, t in enumerate(pool.map(self.interpolate_trajectories, range(self.ntraj))):
                self.z_interpolated[i, :] = t
        else:
            for t in tqdm.tqdm(range(self.ntraj)):
                self.z_interpolated[t, :] = self.trajectories[t, np.digitize(self.time_uniform,
                                                              self.trajectories[t, :, 0], right=False) - 1, 1]
                #self.z_interpolated[t, :] = self.interpolate_trajectories(t, noise=noise)

    def interpolate_trajectories(self, t, noise=0):

        interpolated = np.zeros([self.time_uniform.size])
        for i, x in enumerate(self.time_uniform):
            time_index = np.argmin(np.abs(x - self.trajectories[t, :, 0]))
            if x - self.trajectories[t, time_index, 0] < 0:
                time_index -= 1
            interpolated[i] = self.trajectories[t, time_index, 1]

        # self.pbar.update(1)

        return interpolated + np.random.normal(scale=noise, size=self.time_uniform.size)
        #     self.z_interpolated[t, i] = self.trajectories[t, time_index, 1]
        # self.z_interpolated[t, :] += np.random.normal(scale=noise, size=self.time_uniform.size)

    def calculate_msd(self, ensemble=False):
        """ Calculate mean squared displacement of time series

        :param ensemble: if True, calculate the ensemble msd

        :type ensemble: bool
        """

        print('Calculating MSD...', end='', flush=True)
        start = timer.time()
        self.msd = timeseries.msd(self.z_interpolated.T[..., np.newaxis], 0, ensemble=ensemble, nt=self.nt).T
        print('Done in %.3f seconds' % (timer.time() - start))

    def plot_trajectory(self, n, show=False, save=True, savename='ctrw_trajectory.pdf'):
        """ Plot a CTRW trajectory

        :param n: Trajectory number
        :param show: show plot
        :param save: save plot under savename
        :param savename: name under which to save plot

        :type n: int
        :type show: bool
        :type save: bool
        :type savename: str
        """

        plt.figure()
        # plt.plot(self.trajectory_hops[n, :, 0] / 1000000, self.trajectory_hops[n, :, 1], linewidth=2)
        plt.plot(self.trajectories[n, :, 0] / 1000000, self.trajectories[n, :, 1], linewidth=2)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.xlabel('Time (ms)', fontsize=14)
        plt.ylabel('$z$-coordinate (nm)', fontsize=14)
        plt.tight_layout()

        if show:
            plt.show(block=True)
        if save:
            plt.savefig(savename)

    def bootstrap_msd(self, nboot=200, fit_power_law=False, fit_linear=False):

        if fit_power_law or fit_linear:
            self.fit_parameters = np.zeros([nboot, 2])

        # The average MSD is a collective property, so each bootstrap trial should be an average of self.ntrials
        # ranodomly reconstructed simulated trajectories

        print('Bootstrapping...')

        # This is not any faster than the serial version for some reason
        # if self.nt > 1:
        #     with Pool(self.nt) as pool:
        #         self.bootstraps = np.array(pool.map(self.bootstrap_trial, range(nboot)))
        # else:
        self.bootstraps = np.zeros([nboot, self.msd.shape[1]])
        for i in tqdm.tqdm(range(nboot)):
            self.bootstraps[i, :] = self.bootstrap_trial()

        for i in range(nboot):
            if fit_power_law:
                self.fit_parameters[i, :] = self.fit_power_law(self.bootstraps[i, :])
            if fit_linear:
                self.fit_parameters[i, :] = self.fit_line(self.bootstraps[i, :])

    def bootstrap_trial(self):

        indices = np.random.choice(self.msd.shape[0], size=self.msd.shape[0], replace=True)
        # return self.msd[indices, :].mean(axis=0)  # ~30 % slower than following line (¯\_(ツ)_/¯)

        return self.msd.take(indices, axis=0).mean(axis=0)

    def fit_power_law(self, y, cut=0.6, interactive=True):
        """ Fit power law to MSD curves. Should probably be replaced by MLE
        TODO: weighted fit (need to do error analysis first)

        :param y: y-axis values of MSD curve (x-axis values are values from self.time_uniform
        :param cut: fraction of trajectory to include in fit

        :type y: np.ndarray
        :type cut: float

        :return: Coefficient and exponent in power low of form [coefficient, power]
        """

        self.fit_cut = cut

        start = np.where(y > 1e-6)[0][0]  # find first non-zero value since we will take the log
        end = int(self.fit_cut * len(self.time_uniform))  # fit up until a fraction, cut, of the trajectory
        # print(y[start:end])
        # exit()

        # fit line to linear log plot
        A = Poly_fit.poly_fit(np.log(self.time_uniform[start:end]), np.log(y[start:end]), 1)[-1]

        return [np.exp(A[0]), A[1]]

    def fit_line(self, y, start=0.1, cut=0.5):
        """ Fit line to MSD curve

        :param y: y-axis values of MSD curve (x-axis values are values from self.time_uniform
        :param start: beginning fraction of trajectory to cut out of calculations
        :param cut: fraction of trajectory to include in fit

        :type y: np.ndarray
        :type cut: float

        :return: slope and intercept of fit line [y-intercept, slope]
        """

        self.fit_start = start
        self.fit_cut = cut
        start = int(self.fit_start * len(self.time_uniform))
        end = int(self.fit_cut * len(self.time_uniform))

        A = Poly_fit.poly_fit(self.time_uniform[start:end], y[start:end], 1)[-1]

        return [A[0], A[1]]

    def plot_msd(self, confidence=95, plot_power_law=False, plot_linear=False, show=True, end_frame=None):
        """ Plot averaged mean squared displacement with error bars

        :param confidence: confidence interval for error bars
        :param plot_power_law: if True, fit power law to MSD

        :type confidence: float
        :type plot_power_law: bool
        """

        plt.figure()

        mean = self.msd.mean(axis=0)

        if end_frame is not None:
            mean = mean[:end_frame]
            self.time_uniform = self.time_uniform[:end_frame]

        plt.plot(self.time_uniform, mean, linewidth=2)

        if self.bootstraps is not None:
            error = stats.confidence_interval(self.bootstraps, confidence)
            if end_frame is not None:
                error = error[:, :end_frame]
            plt.fill_between(self.time_uniform, error[1, :] + mean, mean - error[0, :], alpha=0.7)
            self.final_msd = [mean[-1], mean[-1] - error[0, -1], error[1, -1] + mean[-1]]
            print('Estimated MSD: %.2f [%.2f, %.2f]' % (self.final_msd[0], self.final_msd[1], self.final_msd[2]))

        if plot_power_law:
            fit = self.fit_power_law(self.msd.mean(axis=0))
            end = int(self.fit_cut * len(self.time_uniform))
            print('Estimated alpha parameter: %.2f +/- %.2f' % (np.mean(self.fit_parameters[:, 1]),
                                                                np.std(self.fit_parameters[:, 1])))
            plt.plot(self.time_uniform[:end], fitting_functions.power_law(self.time_uniform[:end], fit[0], fit[1]), '--',
                     label='Power law fit')

        if plot_linear:
            fit = self.fit_line(self.msd.mean(axis=0))
            start = int(self.fit_start * len(self.time_uniform))
            end = int(self.fit_cut * len(self.time_uniform))
            print('Estimated slope of line: %.4f +/- %.4f' % (100*np.mean(self.fit_parameters[:, 1]),
                                                                100*np.std(self.fit_parameters[:, 1])))
            plt.plot(self.time_uniform[start:end], fitting_functions.line(fit[1], self.time_uniform[start:end], fit[0]),
                     '--', label='Linear fit')

        if plot_linear | plot_power_law:
            plt.legend()

        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Mean squared displacement (nm$^2$)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()
        #plt.savefig('msd_ctrw.pdf')
        #np.savez_compressed('msd.npz', msd=self.msd.mean(axis=0), error=error, time=self.time_uniform)
        if show:
            plt.show()

    def step_autocorrelation(self):
        """ Calculate autocorrelation of step length and direction
        """

        self.acf = timeseries.step_autocorrelation(self.z_interpolated.T[..., np.newaxis])

    def plot_autocorrelation(self, show=True):
        """ Plot autocorrelation function

        :param show: show plot

        :type show: bool

        :return:
        """

        plt.figure()
        plt.plot(self.time_uniform[:self.acf.shape[1]], self.acf.mean(axis=0))
        plt.xlabel('Time', fontsize=14)
        plt.ylabel('Autocorrelation', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    # ctrw = CTRW(args.steps, args.ntraj, hop_dist=args.hop_length_distribution, dwell_dist=args.dwell_time_distribution,
    #             hop_sigma=args.hop_sigma, alpha=args.alpha, dt=args.dt, nt=args.nthreads, H=args.hurst)
    # ctrw.generate_trajectories(fixed_time=args.fix_time, noise=args.noise, limit=args.upper_limit)
    # ctrw.calculate_msd(ensemble=args.ensemble)
    # ctrw.bootstrap_msd(nboot=args.nboot, fit_power_law=args.fit_power_law, fit_linear=args.fit_line)
    # ctrw.plot_msd(plot_power_law=args.fit_power_law, plot_linear=args.fit_line, show=True)#, end_frame=5000)
    #
    # if args.autocorrelation:
    #     ctrw.step_autocorrelation()
    #     ctrw.plot_autocorrelation(show=True)
    #
    # exit()
    #
    # last = ctrw.msd.mean(axis=0)[int(0.5*ctrw.time_uniform.size)]
    # CI = stats.confidence_interval(ctrw.bootstraps, 95)[:, int(0.5*ctrw.time_uniform.size)]
    # print("MSD at 50 %% of time: %.2f 95 %% CI [%.2f, %.2f]" % (last, last - CI[0], CI[1] + last))
    #
    # # run below with args.ensemble = False, fit_line = True
    # ctrw.calculate_msd(ensemble=True)
    # ctrw.bootstrap_msd(nboot=args.nboot, fit_power_law=True)
    # ctrw.plot_msd(plot_power_law=True, show=False)
    #
    # last = ctrw.msd.mean(axis=0)[-1]
    # CI = stats.confidence_interval(ctrw.bootstraps, 95)[:, -1]
    # print("MSD at 100 %% of time: %.2f 95 %% CI [%.2f, %.2f]" % (last, last - CI[0], CI[1] + last))
    # # plt.show()

    # for ageing demonstration
    steps = [500, 1000, 2000]
    walks = []

    plt.figure()
    for i in steps:
        ctrw = CTRW(i, args.ntraj, hop_dist=args.hop_length_distribution,
                    dwell_dist=args.dwell_time_distribution,
                    hop_sigma=args.hop_sigma, alpha=args.alpha)
        ctrw.generate_trajectories(fixed_time=True)
        ctrw.calculate_msd(ensemble=True)
        plt.plot(ctrw.time_uniform, ctrw.msd.mean(axis=0), linewidth=2, label='Length = %s ns' % i)

    plt.legend()
    plt.xlabel('Time (ns)', fontsize=14)
    plt.ylabel('Mean squared displacement (nm$^2$)', fontsize=14)
    plt.gcf().get_axes()[0].tick_params(labelsize=14)
    plt.tight_layout()
    plt.show()
    exit()

    # for plotting MSDs using a bunch of different dwell time limits
    # limits = [800, 1600, 3200, 6400, 12800, 25600, 51200, 102800]
    # walks = []
    # ctrw = CTRW(args.steps, args.ntraj, hop_dist=args.hop_length_distribution, dwell_dist=args.dwell_time_distribution,
    #             hop_sigma=args.hop_sigma, alpha=args.alpha)
    #
    # plt.figure()
    # for i in limits:
    #     ctrw.generate_trajectories(limit=i)
    #     ctrw.calculate_msd(ensemble=True)
    #     plt.plot(ctrw.time_uniform, ctrw.msd.mean(axis=0), linewidth=2, label='Limit = %s ns' % i)
    #
    # plt.legend()
    # plt.xlabel('Time (ns)', fontsize=14)
    # plt.ylabel('Mean squared displacement (nm$^2$)', fontsize=14)
    # plt.gcf().get_axes()[0].tick_params(labelsize=14)
    # plt.tight_layout()
    # plt.show()
    # exit()

    # for generated MSDs with varied amount of noise (with all else fixed).
    # noise = [1, 0.1, 0]
    # walks = []
    # ctrw = CTRW(args.steps, args.ntraj, hop_dist=args.hop_length_distribution, dwell_dist=args.dwell_time_distribution,
    #             hop_sigma=args.hop_sigma, alpha=args.alpha)
    #
    # plt.figure()
    # for i, n in enumerate(noise):
    #     np.random.seed(1)  # set random seed so trajectory generation will always be the same
    #     ctrw.generate_trajectories(noise=n, nt=args.nthreads)
    #     ctrw.calculate_msd(ensemble=args.ensemble)
    #     #plt.plot(ctrw.time_uniform, ctrw.z_interpolated[1, :] - i*5, label='$\sigma_{noise}$ = %s' % n) # when plotting trajectories
    #     plt.plot(ctrw.time_uniform, ctrw.msd.mean(axis=0), linewidth=2, label='$\sigma_{noise}$ = %s' % n) # when plotting MSD
    #
    # plt.legend()
    # plt.xlabel('Time (ns)', fontsize=14)
    # #plt.ylabel('$z$-position (nm)', fontsize=14)  # when plotting trajectories
    # plt.ylabel('Mean squared displacement (nm$^2$)', fontsize=14)  # when plotting MSD
    # plt.gcf().get_axes()[0].tick_params(labelsize=14)
    # plt.tight_layout()
    # plt.show()
    # exit()
