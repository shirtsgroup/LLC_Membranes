#!/usr/bin/env python

import argparse
import numpy as np
import tqdm
import sys
from LLC_Membranes.llclib import timeseries, fitting_functions
from LLC_Membranes.analysis import Poly_fit
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Simulate fractional brownian motion and calculate its MSD')

    parser.add_argument('-hop', '--hop_length_distribution', default='gaussian', help='Functional form of hop length'
                                                                                      'distribution')
    # if more distributions are included, this will need to be more complicated depending what parameters are needed
    parser.add_argument('-hs', '--hop_sigma', default=1, type=float, help='Standard deviation of gaussian distribution '
                                                                          'used for drawing hop lengths')
    parser.add_argument('-dwell', '--dwell_time_distribution', default='power', help='Functional form of dwell time'
                        'distribution (options: power, exponential')
    parser.add_argument('-steps', '--steps', default=1000, type=int, help='Number of steps to take for each'
                                                                          'independent trajectory')
    parser.add_argument('-ntraj', '--ntraj', default=100, type=int, help='Number of independent ctrw trajectories.')
    parser.add_argument('-ensemble', '--ensemble', action="store_true", help='Calculate MSD as ensemble average')
    parser.add_argument('-power_law', '--fit_power_law', action="store_true", help='Fit MSD to a power law')
    parser.add_argument('-alpha', '--alpha', default=0.5, type=float, help='Anomalous exponent')
    parser.add_argument('-lambda', '--lambda', default=0.5, type=float, help='Exponential decay rate')
    parser.add_argument('-dt', '--dt', default=1, type=float, help='Discrete time step for fixed length simulations')
    parser.add_argument('-fix_time', '--fix_time', action="store_true", help='Fix the total time of simulated '
                        'trajectories. Total length will be steps*dt')

    return parser


def random_exponential_dwell(lam, size=1):
    """ Randomly draw from an exponential distribution
    See Appendix D of https://epubs.siam.org/doi/abs/10.1137/070710111

    :param lam: rate of decay
    :param size: number of random draws to perform

    :type lam: float
    :type size: int

    :return: array of random draws
    """
    return -np.log(1 - np.random.uniform(0, 1, size=size)) / lam


def random_power_law_dwell(alpha, ll=0.1, size=1):
    """ Randomly draw from a power law distribution of form t**-alpha
    See Appendix D of https://epubs.siam.org/doi/abs/10.1137/070710111

    :param alpha: -1 - (anomalous exponent)
    :param ll: lower limit of distribution.
    :param size: number of random draws to perform

    :type: alpha: float
    :type ll: float
    :type size: int

    :return: array of random power law draws
    """

    return ll * (1 - np.random.uniform(0, 1, size=size)) ** (-1 / ((1 + alpha) - 1))


class CTRW(object):

    def __init__(self, length, ntraj, hop_dist='gaussian', dwell_dist='power', hop_sigma=1, alpha=0.5, lamb=0.5,
                 padding=10, dt=1):
        """ Initialize simulation of a continuous time random walk

        :param length: length of each simulated trajectory. If you fix the number of steps, this equals the number of
        steps. If you fix the time, the total length of the simulation is length * dt
        :param ntraj: number of independent trajectories to generate
        :param hop_dist: Name of probability distribution function used to generate random hop lengths
        :param dwell_dist: Name of probability distribution function used to generate random dwell times
        :param hop_sigma: Sigma for Gaussian hop_dist random draws
        :param alpha: Anomalous exponent for power law random draws
        :param lamb: rate of decay for exponential random draws
        :param padding: multiplies number of discrete time points used to interpolate trajectories
        :param dt: time step for fixed time simulations

        :type length: int
        :type ntraj: int
        :type hop_dist: str
        :type dwell_dist: str
        :type hop_sigma: float
        :type alpha: float
        :type lamb: float
        :type padding: int
        :type dt: float
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

        self.trajectories = np.zeros([self.ntraj, self.nsteps, 2])
        self.trajectory_hops = np.zeros([self.ntraj, 2 * self.nsteps - 1, 2])  # for visualization
        self.time_uniform = None
        self.z_interpolated = np.zeros([self.ntraj, self.nsteps*10])  # separate from time_uniform to save memory

        self.msd = None
        self.fit_parameters = None

    def generate_trajectories(self, fixed_time=False):

        if fixed_time:
            self.fixed_time_trajectories()
        else:
            self.fixed_steps_trajectories()

    def fixed_time_trajectories(self):

        length = self.dt * self.nsteps
        print(length)

        for i in tqdm.tqdm(range(self.ntraj)):
            time = [0]

            while np.cumsum(time)[-1] < length:

                if self.hop_distribution == 'gaussian' or self.hop_distribution == 'Gaussian':
                    z_position = np.random.normal(loc=0, scale=self.hop_sigma)
                else:
                    sys.exit('Please enter a valid hop distance probability distribution')

                # hop at random time intervals according to one of the following PDFs
                if self.dwell_distribution == 'exponential':
                    time.append(random_exponential_dwell(self.lamb))
                elif self.dwell_distribution == 'power':
                    time.append(random_power_law_dwell(1 + self.alpha))
                else:
                    sys.exit('Please enter a valid dwell time probability distribution')

            print(length, np.cumsum(time)[-2:])
            exit()

            # make uniform time intervals with the same interval for each simulated trajectory
            max_time = np.min(self.trajectories[:, -1, 0])
            self.time_uniform = np.linspace(0, max_time, self.nsteps * 10)
            for t in tqdm.tqdm(range(self.ntraj)):
                for i, x in enumerate(self.time_uniform):
                    time_index = np.argmin(np.abs(x - self.trajectories[t, :, 0]))
                    if x - self.trajectories[t, time_index, 0] < 0:
                        time_index -= 1
                    self.z_interpolated[t, i] = self.trajectories[t, time_index, 1]

            self.trajectories[i, :, 1] = z_position - z_position[0]  # make initial z equal to 0

    def fixed_steps_trajectories(self):
        """ Generate CTRW trajectories using a fixed number of steps
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
                time = random_exponential_dwell(self.lamb, size=self.nsteps)
            elif self.dwell_distribution == 'power':
                time = random_power_law_dwell(1 + self.alpha, size=self.nsteps)
            else:
                sys.exit('Please enter a valid dwell time probability distribution')

            time = np.cumsum(time)  # accumulate dwell times

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
        for t in tqdm.tqdm(range(self.ntraj)):
            for i, x in enumerate(self.time_uniform):
                time_index = np.argmin(np.abs(x - self.trajectories[t, :, 0]))
                if x - self.trajectories[t, time_index, 0] < 0:
                    time_index -= 1
                self.z_interpolated[t, i] = self.trajectories[t, time_index, 1]

        # plt.plot(self.trajectory_hops[0, :, 0], self.trajectory_hops[0, :, 1])
        # plt.plot(self.time_uniform, self.z_interpolated[0, :])
        # plt.show()
        # exit()

    def calculate_msd(self, ensemble=False):
        """ Calculate mean squared displacement of time series

        :param ensemble: if True, calculate the ensemble msd

        :type ensemble: bool
        """

        print('Calculating MSD...', end='', flush=True)
        self.msd = timeseries.msd(self.z_interpolated.T[..., np.newaxis], 0, ensemble=ensemble).T
        print('Done!')

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

    # def fit_msd(self, start=0, end=-1):
    #     """Interactively fit a line to the mean squared displacement curve
    #
    #     """
    #
    #     self.startfit = start
    #     self.endfit = end
    #     dt = self.time_uniform[1] - self.time_uniform[0]
    #
    #     plt.figure(2)
    #     fit = 0
    #     while fit == 0:
    #
    #         self.yfit = Poly_fit.poly_fit(self.time_uniform[self.startfit:self.endfit],
    #                                       self.msd.mean(axis=0)[self.startfit:self.endfit], 1)[0]
    #
    #         plt.plot(self.time_uniform[self.startfit:self.endfit], self.yfit, '--', color='black',
    #                  label='Linear Fit')
    #         plt.plot(self.time_uniform, self.msd.mean(axis=0), label='MSD')
    #
    #         plt.ylabel('MSD ($nm^2$)', fontsize=14)
    #         plt.xlabel('time (ns)', fontsize=14)
    #         plt.gcf().get_axes()[0].tick_params(labelsize=14)
    #         plt.legend(loc=2)
    #         plt.tight_layout()
    #         plt.ion()
    #         plt.show()
    #         fit = int(input("Type '1' if the fit looks good: "))
    #         if fit != 1:
    #             print('Press enter to following prompts to leave as is')
    #             self.startfit = float(input("Time to start fit (ns): ") or self.startfit)
    #             self.endfit = float(input("Time to stop fit (ns): ") or self.endfit)
    #             self.startfit = int(self.startfit / dt)  # convert time to index in t.time
    #             self.endfit = int(self.endfit / dt)
    #             plt.clf()
    #         else:
    #             plt.close(2)
    #
    # def bootstrap_msd(self, nboot=200):
    #
    #     # The average MSD is a collective property, so each bootstrap trial should be an average of self.ntrials
    #     # ranodomly reconstructed simulated trajectories
    #     self.bootstrapped_msd = np.zeros([nboot, self.msd.shape[1]])
    #     for i in range(nboot):
    #         indices = np.random.choice(self.msd.shape[0], size=self.msd.shape[0], replace=True)
    #         self.bootstrapped_msd[i, :] = self.msd[indices, :].mean(axis=0)
    #
    #     slopes = []
    #     for i in range(nboot):
    #         A = Poly_fit.poly_fit(self.time_uniform[self.startfit:self.endfit],
    #                               self.bootstrapped_msd[i, self.startfit:self.endfit], 1)[-1]
    #         slopes.append(A[1])
    #
    #     self.D = [np.mean(slopes) / (2 * 1 * 10 ** 9),
    #               np.std(slopes) / (2 * 1 * 10 ** 9)]  # divide by dimension and converted to m^2/s
    def fit_power_law(self, cut=0.4, interactive=True):
        """ Fit power law to MSD curves
        TODO: weighted fit (need to do error analysis first)
        :return: Coefficient and exponent in power low of form [coefficient, power]
        """

        end = int(cut * len(self.time_uniform))  # fit up until a fraction, cut, of the trajectory
        print(end)

        # fit line to linear log plot
        A = Poly_fit.poly_fit(np.log(self.time_uniform[1:end]), np.log(self.msd.mean(axis=0)[1:end]), 1)[-1]

        self.fit_parameters = [np.exp(A[0]), A[1]]

    def plot_msd(self, CI=95, nerrorbars=50, fracshow=0.5, save=True, plot_power_law=False):
        """ Plot averaged mean squared displacement with error bars

        :param CI: confidence interval for error bars
        :param nerrorbars: show this many error bars
        :param fracshow: fraction of MSD to plot
        :param save: save the figure and the msd raw data

        :type CI: float
        :type nerrorbars: int
        :type fracshow: float between 0 and 1
        :type save: bool

        """

        plt.figure()
        # error = confidence_interval(self.bootstrapped_msd, CI)
        # plt.errorbar(self.time_uniform, self.msd.mean(axis=0), yerr=error,
        #              errorevery=self.msd.shape[1] // nerrorbars,
        #              linewidth=2, elinewidth=2)
        plt.plot(self.time_uniform, self.msd.mean(axis=0))

        if plot_power_law:
            self.fit_power_law()
            print('Estimated alpha parameter: %.2f' % (self.fit_parameters[1]))
            plt.plot(self.time_uniform, fitting_functions.power_law(self.time_uniform, self.fit_parameters[0],
                                                            self.fit_parameters[1]), '--')

        #plt.title('Diffusivity: %1.2e $\pm$ %1.2e m$^2$/s' % (self.D[0], self.D[1]))
        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Mean squared displacement (nm$^2$)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()
        #plt.savefig('msd_ctrw.pdf')
        #np.savez_compressed('msd.npz', msd=self.msd.mean(axis=0), error=error, time=self.time_uniform)
        plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    ctrw = CTRW(args.steps, args.ntraj, hop_dist=args.hop_length_distribution, dwell_dist=args.dwell_time_distribution,
                hop_sigma=args.hop_sigma, alpha=args.alpha, dt=args.dt)
    ctrw.generate_trajectories(fixed_time=args.fix_time)
    ctrw.calculate_msd(ensemble=args.ensemble)
    ctrw.plot_msd()

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