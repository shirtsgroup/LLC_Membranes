#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.analysis import Poly_fit
from LLC_Membranes.llclib import timeseries, fitting_functions
import tqdm
import fbm


def initialize():

    parser = argparse.ArgumentParser(description='Simulate fractional brownian motion and calculate its MSD')

    parser.add_argument('-H', '--hurst', default=0.5, type=float, help='Hurst exponent for fbm simulation. 2H=alpha '
                        'where alpha is the anomalous exponent. H=0.5 is brownian motion')
    parser.add_argument('-fsteps', '--fbmsteps', default=1000, type=int, help='Number of steps to take along fractional'
                                                                              'Brownian motion trajectory')
    parser.add_argument('-nftraj', '--nfbmtraj', default=100, type=int, help='Number of independent fractional brownian'
                                                                             ' motion trajectories.')
    parser.add_argument('-fdt', '--ftimestep', default=1, type=float, help='Time step between steps in fractional '
                        'brownian motion simulated trajectory (arbitrary units).')
    parser.add_argument('-acf', '--autocorrelation', action="store_true", help='Plot step autocorrelation function')
    parser.add_argument('-acov', '--autocovariance', action="store_true", help='Plot autocovariance of fractional'
                                                                               'Gaussian noise')
    parser.add_argument('-ensemble', '--ensemble', action="store_true", help='Calculate MSD as ensemble average')
    parser.add_argument('-power_law', '--power_law', action="store_true", help='Fit MSD to a power law')
    parser.add_argument('-scale', '--scale', default=None, type=float, help='Sigma of hop length distribution')

    return parser


class FractionalBrownianMotion(object):

    def __init__(self, nsteps, hurst, ntraj=1, length=1, method="daviesharte", scale=1, progress=True):
        """ Generate trajectories exhibiting fractional brownian motion

        :param nsteps: number of steps in trajectory
        :param hurst: hurst exponent (2H=alpha where alpha is anomalous exponent).
        :param ntraj: number of independent trajetories to generate
        :param length: total length of simulation trajectory
        :param method: method of fractional brownian motion simulation
        :param scale: standard deviation in step size. A bigger scale means bigger steps.
        """

        self.hurst = hurst
        self.ntraj = ntraj
        self.nsteps = nsteps

        self.fbm = fbm.FBM(self.nsteps, self.hurst, length=length, method=method)
        self.time = fbm.times(self.nsteps, length=length)

        self.trajectories = np.zeros([nsteps + 1, ntraj])

        if progress:
            print('Generating FBM trajectories...')
        for i in tqdm.tqdm(range(ntraj), disable=(not progress)):
            self.trajectories[:, i] = self.fbm.fbm()*scale

        self.msd = None
        self.fit_parameters = None
        self.acf = None
        self.acov = None

    def plot_fbm(self):
        """ Plot a fractional brownian motion trajectory using parameters in __init__
        """

        plt.figure()
        plt.plot(self.time, self.trajectories[:, 0], linewidth=2)
        plt.xlabel('Time', fontsize=14)
        plt.ylabel('Position', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        plt.show()

    def calculate_msd(self, ensemble=False):
        """ Calculate mean squared displacement of time series

        :param ensemble: if True, calculate the ensemble msd

        :type ensemble: bool
        """

        print('Calculating MSD...', end='', flush=True)
        self.msd = timeseries.msd(self.trajectories[:, :, np.newaxis], 0, ensemble=ensemble)
        print('Done!')

    def fit_power_law(self):
        """ Fit power law to MSD curves
        TODO: weighted fit (need to do error analysis first)
        :return: Coefficient and exponent in power low of form [coefficient, power]
        """

        # fit line to linear log plot
        A = Poly_fit.poly_fit(np.log(self.time[1:]), np.log(self.msd.mean(axis=1)[1:]), 1)[-1]

        self.fit_parameters = [np.exp(A[0]), A[1]]

    def plot_msd(self, plot_power_law=True, show=True):
        """ Plot mean squared displacement with optional power law fit

        :param plot_power_law: plot power law fit to data
        :param show: show plot

        :type plot_power_law: bool
        :type show: bool

        :return:
        """

        plt.figure()
        plt.plot(self.time, self.msd.mean(axis=1), linewidth=2)

        if plot_power_law:
            self.fit_power_law()
            print('Estimated Hurst parameter: %.2f' % (self.fit_parameters[1] / 2))
            plt.plot(self.time, fitting_functions.power_law(self.time, self.fit_parameters[0],
                                                            self.fit_parameters[1]), '--')

        plt.xlabel('Time', fontsize=14)
        plt.ylabel('MSD', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()

    def step_autocorrelation(self):
        """ Calculate autocorrelation of step length and direction
        """

        self.acf = timeseries.step_autocorrelation(self.trajectories[..., np.newaxis])

    def plot_autocorrelation(self, show=True):
        """ Plot autocorrelation function

        :param show: show plot

        :type show: bool

        :return:
        """

        plt.figure()
        plt.plot(self.time[:self.acf.shape[1]], self.acf.mean(axis=0))
        plt.xlabel('Time', fontsize=14)
        plt.ylabel('Autocorrelation', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()

    def step_autocovariance(self):
        """ Calculate autocovariance of fractional gaussian noise in the trajectories (i.e. the step lengths)
        """

        self.acov = timeseries.autocov((self.trajectories[1:, :] - self.trajectories[:-1, :]).T)

    def analytical_autocovariance_fbm(self):

        k = np.arange(len(self.time))

        return 0.5 * (np.abs(k - 1) ** (2 * self.hurst) - 2 * k ** (2 * self.hurst) + (k + 1) ** (2 * self.hurst))

    def plot_autocovariance(self, show=False, analytical=True):
        """ Plot autocovariance function

        :param show: show plot
        :param analytical: plot analytical autocovariance function for fractional gaussian noise

        :type show: bool
        :type analytical: bool

        """

        plt.figure()

        plt.plot(self.time[:-1], self.acov, label='Simulated', color='blue', linewidth=3)

        if analytical:
            plt.plot(self.time, self.analytical_autocovariance_fbm(), label='Analytical', linewidth=1, color='black')
            plt.legend()

        plt.xlabel('Time Lag', fontsize=14)
        plt.ylabel('Autocovariance', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    fbm = FractionalBrownianMotion(args.fbmsteps, args.hurst, ntraj=args.nfbmtraj,
                                   length=(args.ftimestep*args.fbmsteps), method="daviesharte", scale=args.scale)

    fbm.calculate_msd(ensemble=args.ensemble)

    if args.autocorrelation:
        fbm.step_autocorrelation()
        print(fbm.acf[:, 1].mean())
        fbm.plot_autocorrelation(show=False)

    if args.autocovariance:
        fbm.step_autocovariance()
        print(fbm.acov[1])
        fbm.plot_autocovariance(show=False, analytical=False)
    exit()

    fbm.plot_msd(plot_power_law=args.power_law, show=True)
