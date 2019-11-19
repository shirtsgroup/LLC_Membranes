#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import levy
import argparse
import tqdm


def initialize():

    parser = argparse.ArgumentParser(description='Determine the mean first passage time for a particle performing a'
                                                 'random walk.')

    parser.add_argument('-n', '--nwalks', default=1000, type=int, help='Number of random walk trajectories to generate')
    parser.add_argument('-L', '--length', default=10, type=float, help='Length of domain in which particle travels')
    parser.add_argument('-x0', '--x0', default=5, type=float, help='Initial particle position within domain')
    parser.add_argument('-sigma', '--sigma', default=0.5, type=float, help='Standard deviation of hop length per unit'
                        'time. This effectivley determines the diffusion constant.')
    parser.add_argument('-dt', '--time_step', default=0.01, type=float, help='Time between draws from hop distribution.'
                        ' The smaller, the better.')

    # plotting
    parser.add_argument('-maxt', '--maxt', default=100, type=float, help='Maximum passage time to include in plots.')
    parser.add_argument('-savename', '--savename', help='Name of saved plots. cdf and pdf will be added to end of '
                                                        'names of saved plots.')

    return parser


class CDF:

    def __init__(self, data):
        """ Generate an emperical cumulative distribution function from data

        :param data: x-values of data in no particular order

        :type data: list
        """

        self.xs = np.array(sorted(data))
        self.N = float(len(self.xs))
        self.ys = np.arange(1, self.N + 1) / self.N

    def cdf(self, x):
        """ Callable cumulative emperical distribution function
        :param x: array of x-values at which to evaluate cumulative emperical distribution function
        :return:
        """

        if type(x) is np.float64:
            x = np.array([x])

        ndx = [np.argmin(np.abs(self.xs - x[i])) for i in range(x.size)]

        return self.ys[ndx]


class MFPT:

    def __init__(self, L, x0):
        """ Release a particle at x0 and let it perform a random walk until it is absorbed at 0 or L

        :param L: length
        :param x0: initial particle position

        :type L: float
        :type x0: float
        """

        # system geometry and release point
        self.L = L
        self.x0 = x0

        # Initializing variables for usage later
        self.passage_times = None
        self.dt = None
        self.sigma = None
        self.trajectories = []

    def simulate_passage_time(self, nwalks, sigma, dt, nplot_random=0, nt=1):
        """
        :param nwalks: number of trajectories to simulate
        :param sigma: width of hop distribution. Related to diffusion constant.
        :param dt: time step
        :param nplot_random: number of random trajectories to plot. If 0, a plot won't be made

        :type nwalks: int
        :type sigma: float
        :type dt: float
        :type nplot_random: int
        """

        self.passage_times = np.zeros(nwalks)
        self.dt = dt
        self.sigma = sigma
        jump_sigma = np.sqrt(dt) * sigma

        for w in tqdm.tqdm(range(self.passage_times.size), unit='walks'):

            walk = [self.x0]
            while 0 < walk[-1] < self.L:

                walk.append(walk[-1] + jump_sigma * levy.random(2, 0, mu=0))  # faster than scipy

            self.passage_times[w] = (len(walk) - 1) * self.dt  # subtract 1 so time zero isn't counted

            self.trajectories.append(walk)

        if nplot_random > 0:

            # plot random trajectories
            ndx = np.random.choice(nwalks, size=nplot_random)
            for i in ndx:
                plt.plot(self.dt*np.arange(int(self.passage_times[i] / self.dt) + 1), self.trajectories[i])

            plt.show()

    def plot_cdf(self, max_t=200, show=True, savename=False):
        """ Plot the analtyical versus empirical cumulative distribution function of first passage times

        :param max_t: Largest value of time to show
        :param show: show the plot once it's made
        :param savename: if not None, save the figure by this name

        :type max_t: float
        :type show: bool
        :type savename: NoneType or str
        """

        plt.figure()
        time = self.dt*np.arange(1, int(max_t / self.dt)).astype(float)

        plt.plot(time, self._analytical_passage_time_cdf(time), lw=2, label='analytical')
        plt.plot(time, self._empirical_cdf(time), lw=2, label='empirical')

        # formatting
        plt.legend(fontsize=14)
        plt.xlabel('Time', fontsize=14)
        plt.ylabel('Cumulative Density', fontsize=14)
        plt.tick_params(labelsize=14)
        plt.tight_layout()

        if savename is not None:
            plt.savefig('%s_cdf.pdf' % savename)

        if show:
            plt.show()

    def plot_passage_time_distribution(self, bins=50, max_t=100, show=True, savename=None):
        """ Plot the distribution of first passage times against the analytical PDF

        :param bins: number of bins to use to discretize empirical data
        :param max_t: maxmimum passage time to plot
        :param show: show the plot when done
        :param savename: if not None, save the figure by this name

        :type bins: int
        :type max_t: float
        :type show: bool
        :type savename: NoneType or str
        """

        plt.figure()  # create new figure

        # get PDF by taking derivative of CDF
        time = self.dt * np.arange(1, int(max_t / self.dt)).astype(float)
        cdf = self._analytical_passage_time_cdf(time)
        time = time[cdf > 0]  # there is a giant drop near zero
        cdf = cdf[cdf > 0]
        dx = time[1] - time[0]
        deriv = np.diff(cdf) / dx
        plt.plot(time[:-1], deriv, '--', label='Approximate Analytical PDF', lw=2, color='black')

        # histogram empirical measurements
        heights, edges = np.histogram(self.passage_times, bins=bins, range=(0, max_t), density=True)
        heights = heights / cdf[-1]  # normalize by area under theoretical curve up to max_t
        bin_width = edges[1] - edges[0]
        bin_centers = [i + bin_width/2 for i in edges[:-1]]

        plt.bar(bin_centers, heights, bin_width, label='Observations')

        print('Empirical mean first passage time: %.2f' % self.passage_times.mean())
        print('Approximate analytical mean first passage time: %.2f' % np.average(time[:-1], weights=deriv))

        plt.legend(fontsize=14)
        plt.xlabel('Time', fontsize=14)
        plt.ylabel('Probability', fontsize=14)
        plt.tick_params(labelsize=14)
        plt.tight_layout()

        if savename is not None:
            plt.savefig('%s_pdf.pdf' % savename)

        if show:
            plt.show()

    def _analytical_passage_time_cdf(self, time, nterms=100):
        """ Equation 6 from : https://journals.aps.org/pre/pdf/10.1103/PhysRevE.73.046104
        """

        f = np.zeros_like(time)
        for i, t in enumerate(time):
            summation = 0
            for j in range(1, nterms + 1, 2):  # only odd terms are nonzero
                summation += self._analytical_cdf_term(j, t)
            f[i] = 1 - (2 / np.pi) * summation

        return f

    def _analytical_cdf_term(self, j, t):
        """ Terms from summation in CDF

        :param j:
        :param t:
        :return:
        """

        try:

            return ((1 - np.cos(j * np.pi)) / j) * np.sin(j * np.pi * self.x0 / self.L) * \
                    np.exp(-(j * np.pi * self.sigma / self.L) ** 2 * t)

        except FloatingPointError:  # at large t, the exponential term gets REALLY small causing an underflow
            return 0

    def _empirical_cdf(self, x):
        """ Use CDF class to create an empirical cumulative distribution function (ECDF) from the observed first
        passage times

        :param x: discrete values at which to get approximate ECDF values

        :type x: list or numpy.ndarray
        """

        return CDF(self.passage_times.tolist()).cdf(x)


if __name__ == "__main__":

    args = initialize().parse_args()

    random_walks = MFPT(args.length, args.x0)
    random_walks.simulate_passage_time(args.nwalks, args.sigma, args.time_step, nplot_random=0)
    random_walks.plot_cdf(100, show=False, savename=args.savename)
    random_walks.plot_passage_time_distribution(50, max_t=100, show=True, savename=args.savename)
