#!/usr/bin/env python

""" Generate Linear Fractional Stable Noise
"""

import numpy as np
from scipy.stats import levy_stable
import sys
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import timeseries, stats
import tqdm


class FLM:

    def __init__(self, H, alpha, m=256, M=6000, C=1, N=2**12, n=1, scale=1):
        """ Generate realizations of fractional levy motion, also know as linear fractional stable motion

        :param H: Hurst parameter. Also known as the self-similarity parameter
        :param alpha: the tail-exponent of the stable distribution (between 0 and 2). Lower alpha = heavier tails
        :param m: 1/m is the mesh size
        :param M: kernel cut-off parameter
        :param C: scale parameter
        :param N: size of sample
        :param scale:
        """

        self.H = H
        self.alpha = alpha
        self.m = m
        self.M = M
        self.N = N
        self.scale = C
        self.Na = m * (M + N)

        if alpha < 0 or alpha > 2:
            sys.exit('Alpha must be greater than 0 and less than or equal to 2!')

        mh = 1 / m
        d = H - 1 / self.alpha
        t0 = np.linspace(mh, 1, m) ** d
        t1 = np.linspace(1 + mh, M, int((M - (1 + mh)) / mh) + 1)
        t1 = t1 ** d - (t1 - 1) ** d
        self.A = mh ** (1 / alpha) * np.concatenate((t0, t1))
        self.scale *= (np.abs(self.A) ** alpha).sum() ** (-1 / alpha)
        self.A *= self.scale
        self.A = np.fft.fft(self.A, n=self.Na)

        self.realizations = None

    def generate_realizations(self, n, truncate=None):
        """ Generate realization of fractional levy motion

        :param n: Number of realizations to generate
        :param truncate: largest allowable fluctuation

        :type n: int
        :type truncate: float or None
        """

        self.realizations = np.zeros([n, self.N])

        for i in tqdm.tqdm(range(n)):

            if self.alpha == 2:
                z = np.random.normal(0, scale=self.scale, size=self.Na)
            else:
                if truncate is None:
                    z = levy_stable.rvs(self.alpha, 0, loc=0, scale=self.scale, size=self.Na)
                else:
                    z = np.zeros(self.Na)
                    # TODO: generate all of them then replace all values with new rvs until all under truncate.
                    for j in tqdm.tqdm(range(self.Na)):
                        rv = levy_stable.rvs(self.alpha, 0, loc=0, scale=self.scale)
                        while rv > truncate:
                            rv = levy_stable.rvs(self.alpha, 0, loc=0, scale=self.scale)
                        z[j] = rv

            z = np.fft.fft(z, self.Na)
            w = np.fft.ifft(z * self.A, self.Na).real

            self.realizations[i, :] = w[:self.N*self.m:self.m]

        print(np.std(self.realizations.flatten()))

    def plot_marginal(self, bounds=(-4, 4), bins=50, show=False):
        """ Plot a histogram of the marginal distribution of increments, with the expect PDF overlayed on top of it

        :param bounds: largest increments to be included in histogram
        :param bins: number of bins in histogram
        :param show: show the plot when done

        :type bounds: tuple of floats
        :type bins: int
        :type show: bool
        """

        x = np.linspace(bounds[0], bounds[1], 1000)

        hist, bin_edges = np.histogram(self.realizations.flatten(), bins=bins, range=bounds, density=True)

        # account for part of PDF that is chopped off. Using density=True makes hist sum to 1
        area_covered = levy_stable.cdf(bounds[1], self.alpha, 0, loc=0, scale=self.scale) - \
                       levy_stable.cdf(bounds[0], self.alpha, 0, loc=0, scale=self.scale)
        hist *= area_covered

        # plot bars. Can't use plt.hist since I needed to modify the bin heights
        bin_width = bin_edges[1] - bin_edges[0]
        bin_centers = [i + bin_width / 2 for i in bin_edges[:-1]]
        plt.figure()
        plt.bar(bin_centers, hist, width=bin_width)
        plt.plot(x, levy_stable.pdf(x, self.alpha, 0, loc=0, scale=self.scale), '--', color='black', lw=2)

        # formatting
        plt.xlabel('Step Size', fontsize=14)
        plt.ylabel('Frequency', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()

    def plot_autocorrelation(self, max_k=25, nboot=200, confidence=68.27, show=False):
        """ Plot autocorrelation function of increments

        :param max_k: maximum lag time to plot
        :param nboot: number of bootstrap trials
        :param confidence: confidence interval of shaded error region (percent)
        :param show: show the plot when done

        :type max_k: int
        :type nboot: int
        :type confidence: float
        :type show: bool
        """

        # calculate acf of each trajectory
        ntraj = self.realizations.shape[0]
        acf = np.zeros([ntraj, self.N - 1])
        for i in range(ntraj):
            acf[i, :] = timeseries.acf(self.realizations[i, :])

        # bootstrap
        boot = np.zeros([nboot, self.N - 1])
        for i in range(nboot):
            ndx = np.random.randint(ntraj, size=ntraj)
            boot[i, :] = acf[ndx, :].mean(axis=0)

        errorbars = stats.confidence_interval(boot, confidence)

        plt.figure()
        avg = boot.mean(axis=0)
        plt.plot(np.arange(self.N - 1), avg, lw=2)
        plt.fill_between(np.arange(self.N - 1), avg + errorbars[0, :], avg - errorbars[1, :], alpha=0.25)

        # formatting
        plt.xlim(-0.5, max_k)
        plt.ylim(-0.6, 1)
        plt.xlabel('Lag Time (steps)', fontsize=14)
        plt.ylabel('Correlation', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()

    def plot_trajectory(self, traj_no, show=False):
        """ Plot trajectory(ies)

        :param traj_no: trajectory number or list of trajetory numbers to plot
        :param show: show the plot when done

        :type traj_no: int or list of ints
        :type show: bool
        """

        if type(traj_no) is int:
            traj_no = [traj_no]

        plt.figure()
        for i in traj_no:
            plt.plot(np.cumsum(self.realizations[i, :]), lw=2)

        plt.xlabel('Time', fontsize=14)
        plt.ylabel('Position', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()

    def plot_msd(self, frac=0.4, nboot=200, confidence=68, show=False):
        """ Calculate and plot the MSD

        :param frac: fraction of MSD plot to show
        :param nboot: number of bootstrap trials
        :param confidence: percent confidence interval
        :param show: show the plot when done

        :type frac: float
        :type nboot: int
        :type confidence: float
        :type show: bool
        """

        plt.figure()

        traj = np.cumsum(self.realizations, axis=1)
        msds = timeseries.msd(traj.T[..., np.newaxis], 0)

        errorbars = timeseries.bootstrap_msd(msds, nboot, confidence=confidence)
        end = int(frac*self.N)
        avg = msds.mean(axis=1)

        plt.fill_between(np.arange(end), avg[:end] + errorbars[0, :end], avg[:end] - errorbars[1, :end], alpha=0.25)
        plt.plot(np.arange(end), avg[:end])
        plt.xlabel('Time', fontsize=14)
        plt.ylabel('MSD (nm$^2$)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()


if __name__ == "__main__":

    np.random.seed(1)
    flm = FLM(0.4, 1.5, C=1)
    flm.generate_realizations(24, truncate=5)
    flm.plot_marginal()
    flm.plot_autocorrelation()
    flm.plot_trajectory([0, 1])
    flm.plot_msd(show=True)

