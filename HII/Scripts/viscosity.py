#!/usr/bin/env python

import mdtraj as md
import argparse
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
from scipy.optimize import curve_fit
from scipy.integrate import cumtrapz
from scipy import stats
import tqdm
import mmap


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Trajectory file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-e', '--edr', default='wiggle.edr', type=str, help='Gromacs portable energy file')
    parser.add_argument('-T', '--temperature', default=300, type=float, help='Temperature of simulation')
    parser.add_argument('-n', '--nsub', default=5, type=int, help='Number of subtrajectories to break simulation'
                        'into. Viscosity is calculated for each subtrajectory and used to generate statistics')
    parser.add_argument('-l', '--load', default=False, help='True or 1 will load everything up until the trajectory'
                        'is broken into sub-intervals. Set this to 2 to load previously calculated viscosities')

    args = parser.parse_args()

    return args


def largest_prime_factor(n):
    i = 2
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
    return n


def prime_factors(n):
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors


def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines


def autocorrelation(x, largest_prime=500):
    """ FFT based autocorrelation function, which is faster than numpy.correlate. Efficiency is key in order to avoid
    headaches.
    :param x : multidimensional numpy array of y values of an equispaced timeseries. [npoints, ntrajectories]
    :param largest_prime : the largest prime factor of array length allowed. The smaller the faster. 1.6M points takes
    about 5 seconds with largest_prime=1000. Just be aware that you are losing data by truncating. But 5-6 data points
    isn't a big deal usually.
    """

    # x -= np.mean(x, axis=0)  # subtract the mean (alternate, but not equivalent way of calculating autocorrelation)

    l = 2*x.shape[0] - 1

    while largest_prime_factor(l) >= largest_prime or l % 2 == 0:
        l -= 1

    x = x[:(l + 1) // 2, :]
    length = x.shape[0]*2 - 1

    fftx = np.fft.fft(x, n=length, axis=0)
    ret = np.fft.ifft(fftx * np.conjugate(fftx), axis=0)
    ret = np.fft.fftshift(ret)

    # divide each point in the autocorrelation function by the number of counts
    auto = ret[length // 2:].real
    auto = np.mean(auto, axis=1)

    n = auto.shape[0]

    # return auto / (x.var()*np.arange(n, 0, -1))  # if subtracting mean, this is right thing to return

    return auto / np.arange(n, 0, -1)


def estimated_autocorrelation(x):
    """ Another way to calculate autocorrelation function. Produces same results as autocorrelation(). This is about
    5 times slower though
    """
    n = len(x)
    variance = x.var()
    x = x - x.mean()
    r = np.correlate(x, x, mode='full')[-n:]
    result = r / (variance * (np.arange(n, 0, -1)))

    return result


def viscosity_fit(t, A, a, tau1, tau2):
    """
    :param p: fit parameters in a list [A, alpha, tau1, tau2]
    :param t: x values (time)
    :return: y values at each x for the function with parameters, p
    """

    return A*a*tau1*(1 - np.exp(-t/tau1)) + A*(1 - a)*tau2*(1 - np.exp(-t/tau2))


class Viscosity(object):
    """
    Calculate the viscosity of system
    """

    def __init__(self, edr, gro, max_lag=10, load=False):
        """
        :param edr: name of energy file (.edr)
        :param gro: name of coordinate file (.gro) (Assumes that system was run nve so volume is constant)
        :param max_lag : maximum time lag to for calculation of autocorrelation function (ps)
        """

        self.kb = 1.38064852 * 10 ** -23  # boltzmann constant [=] m**2 * kg * s**-2 * K**- 1
        self.ksi = 2.837297  # unitless.
        t = md.load(gro)
        box = t.unitcell_vectors[0, :, :]
        self.box_length = box[0, 0]  # all sides the same
        self.volume = np.dot(box[2, :], np.cross(box[0, :], box[1, :]))  # volume of unit cell [=] nm3
        self.volume *= 1 * 10 ** -27  # convert to m3
        self.max_lag = max_lag
        self.autocorr = []
        self.runningintegral = []
        self.fit_params = None
        self.viscosity = []
        self.average_viscosity = 0
        self.std_viscosity = 0
        self.load = load
        self.correction = 0

        if not self.load:

            # get all off diagonal elements of the pressure tensor vs. time using gmx energy
            print('Extracting off-diagonal elements of the pressure tensor and Temperature from %s using gmx_energy ...'
                  % edr, end=' ', flush=True)
            ps = subprocess.Popen(('echo', 'Pres-XY', '\n', 'Pres-XZ', '\n', 'Pres-YX', '\n', 'Pres-YZ', '\n', 'Pres-ZX',
                                   '\n', 'Pres-ZY', '\n', 'Temperature', '\n', 'Pressure', '\n'), stdout=subprocess.PIPE)
            ps.wait()
            p = subprocess.Popen(('gmx', 'energy', '-f', '%s' % edr), stdin=ps.stdout,
                                 stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
            p.wait()

            print('Reading energy.xvg')
            with open('energy.xvg', 'r') as f:
                xvg = []
                for line in tqdm.tqdm(f, total=get_num_lines('energy.xvg')):
                    xvg.append(line)

            # find where the data starts by looking for the line that contains "Pres-ZY" since that was the last option
            # passed to gmx energy
            data_start = 0
            while xvg[data_start].count('Pres-ZY') == 0:
                data_start += 1
            data_start += 1

            self.P = np.zeros([len(xvg) - data_start, 6])
            self.T = np.zeros([len(xvg) - data_start])
            self.total_pressure = np.zeros([len(xvg) - data_start])
            self.time = np.zeros([len(xvg) - data_start])
            for i in range(data_start, len(xvg)):
                data = xvg[i].split()
                self.time[i - data_start] = float(data[0])
                self.P[i - data_start, :] = [float(d) for d in data[3:]]
                self.T[i - data_start] = float(data[1])
                self.total_pressure[i - data_start] = float(data[2])

            self.Tavg = np.mean(self.T)  # [=] K
            self.Pavg = np.mean(self.total_pressure)
            self.dt = self.time[1] - self.time[0]  # [=] ps
            self.P *= 100000  # convert from bar to pascal
            self.P_full = np.copy(self.P)  # save a copy before this gets modified. This will be saved later
            self.prefactor = self.volume / (self.kb * self.Tavg)  # [=] m * s**2 * kg**-1 [=] Pa**-1

            self.nintervals = 0
            while self.time[self.nintervals] < self.max_lag:
                self.nintervals += 1

        else:

            self.zip = np.load('viscosity.npz')
            self.nintervals = self.zip['nintervals']
            self.P = self.zip['P']
            self.prefactor = self.zip['prefactor']
            self.time = self.zip['time']
            self.dt = self.zip['dt']
            self.T = self.zip['T']
            self.total_pressure = self.zip['total_P']

    def break_apart_trajectory(self, nsub):
        """
        Break the trajectory into n equal length trajectories
        :param n: number of subtrajectories to break trajectory into
        """

        if self.load != '2':
            total_length = self.P.shape[0]  # total number of frames
            length_sub = total_length // nsub  # number of frames per subtrajectory
            self.P = self.P[:int(nsub*length_sub), :]  # get rid of extra frames (if any)
            self.P = np.reshape(self.P, (nsub, length_sub, self.P.shape[1]))
            print('Trajectory broken into %d sub-trajectories of length %.1f ps' % (nsub, self.dt*length_sub))

    def calculate(self, nboot=1000):
        """
        Calculate viscosity
        """

        nsub = self.P.shape[0]
        bounds = ([-np.inf, 0, 0, 0], [np.inf, 1, np.inf, np.inf])  # bounds on parameters that will be fit

        if args.load != '2':
            for i in range(nsub):

                print('Calculating viscosity for trajectory # %d' % (i + 1), end='\r', flush=True)
                self.autocorr.append(autocorrelation(self.P[i, ...]).real)

                # calculate running integral of autocorrelation function
                ps_to_s = 1*10**-12  # conversion of picoseconds to seconds
                self.runningintegral.append(ps_to_s * cumtrapz(self.autocorr[i][:self.nintervals], dx=self.dt, initial=0)
                                            * self.prefactor)

                # fit double exponential function to running integral
                solp, pcov = curve_fit(viscosity_fit, self.time[:self.nintervals],
                                       self.runningintegral[i][:self.nintervals], bounds=bounds)

                self.viscosity.append(viscosity_fit(np.inf, solp[0], solp[1], solp[2], solp[3]))

            self.viscosity = np.array(self.viscosity)

            noutliers = self.outliers(0.01)
            if noutliers > 0:
                print('\n%s outlier(s) discarded (alpha = 0.01)' % noutliers)

            self.autocorr = np.mean(np.array(self.autocorr), axis=0)  # The mean autocorrelation function will be plotted
            self.runningintegral = np.mean(np.array(self.runningintegral), axis=0)  # plot mean running integral

            solp, pcov = curve_fit(viscosity_fit, self.time[:self.nintervals],
                                   self.runningintegral[:self.nintervals], bounds=bounds)
            self.fit_params = solp
            self.average_viscosity, self.std_viscosity = self.bootstrap(nboot)
            self.correction = np.mean(self.T) * self.kb * self.ksi / (6 * np.pi * self.average_viscosity * self.box_length)
            self.correction *= 10**9  # convert box_length from nm to m
            self.correction_error = self.correction * (self.std_viscosity / self.average_viscosity)  # propagate viscosity error

        else:
            self.viscosity = self.zip['viscosity']
            self.autocorr = self.zip['autocorr']
            self.fit_params = self.zip['fit_params']
            self.runningintegral = self.zip['running_integral']
            self.average_viscosity = self.zip['average_viscosity']  # outliers should have already been taken care of
            self.std_viscosity = self.zip['std_viscosity']
            self.correction = self.zip['correction']
            self.correction_error = self.zip['correction_error']

        # multiply by 1000 to convert from Pa-s to cP

        print('\nAverage Temperature: %.1f +/- %.1f K' % (np.mean(self.T), np.std(self.T)))
        print('Average Pressure: %.2f +/- %.2f bar' % (np.mean(self.total_pressure), np.std(self.total_pressure)))
        print('Average Viscosity : %.3g +/- %.3g cP' % (1000*self.average_viscosity, 1000*self.std_viscosity))
        print('Correction to Diffusivity = %.2e m^2/s +/- %.2e m^2/s' % (self.correction, self.correction_error))

    def outliers(self, alpha):
        """
        Check for outliers of viscosity calculation using Grubbs' test
        Steps:
        (1) Calculate critical t-statistic
        https://stackoverflow.com/questions/19339305/python-function-to-get-the-t-statistic?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
        (2) Calculate critical G
        https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h1.htm
        :param alpha : probability that point is falsely rejected
        :return:
        """

        outlier = True  # hypothesize that there is an outlier
        noutliers = 0

        while outlier:
            N = len(self.viscosity)  # number of samples
            t = stats.t.ppf(1 - ((alpha/2)/(2*N)), N - 2)
            Gcrit = np.sqrt(t**2 / (N - 2 + t**2))
            Gcrit *= (N - 1)/(N**0.5)
            G = np.abs((self.viscosity - self.viscosity.mean()) / self.viscosity.std())
            potential_outlier = np.amax(G)
            if potential_outlier > Gcrit:
                self.viscosity = np.delete(self.viscosity, np.argmax(G))
                del self.autocorr[np.argmax(G)]
                del self.runningintegral[np.argmax(G)]
                noutliers += 1
            else:
                outlier = False

        return noutliers

    def bootstrap(self, n):

        ndx = np.random.randint(0, len(self.viscosity), size=n)
        trials = [self.viscosity[ndx[i]] for i in range(n)]

        return np.mean(trials), np.std(trials)

    def plot_all(self, save=False):
        """
        Plot all relevant data
        :return:
        """
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(self.time[:10*self.nintervals], self.autocorr[:10*self.nintervals]/np.amax(self.autocorr[:10*self.nintervals]))
        left, bottom, width, height = [0.45, 0.4, 0.42, 0.42]
        ax2 = fig.add_axes([left, bottom, width, height])
        ax2.plot(self.time[:self.nintervals], self.autocorr[:self.nintervals]/np.amax(self.autocorr[:self.nintervals]))
        ax.set_xlabel('Time (ps)', fontsize=14)
        ax.set_ylabel('Autocorrelation', fontsize=14)
        ax2.set_ylim(-0.01, 0.025)
        ax2.set_xlim(0, 10)
        ax.set_xscale('log')
        # plt.tight_layout() # not compatible with inset plot
        if save:
            plt.savefig('autocorrelation_function.py')

        plt.figure()
        plt.plot(self.time[:self.nintervals], viscosity_fit(np.array(self.time[:self.nintervals]), self.fit_params[0],
                 self.fit_params[1], self.fit_params[2], self.fit_params[3]), '--', label='Double Exponential Fit')
        plt.plot(self.time[:self.nintervals], self.runningintegral[:self.nintervals], label='Running Integral')
        plt.xlabel('Time (ps)', fontsize=14)
        plt.ylabel('$\eta(t)$')
        plt.legend()
        plt.tight_layout()
        if save:
            plt.savefig('running_integral.png')

        plt.show()

    def save(self):

        np.savez_compressed('viscosity.npz', nintervals=self.nintervals, P=self.P_full, prefactor=self.prefactor,
                            time=self.time, dt=self.dt, viscosity=self.viscosity, autocorr=self.autocorr,
                            fit_params=self.fit_params, running_integral=self.runningintegral,
                            average_viscosity=self.average_viscosity, std_viscosity=self.std_viscosity,
                            correction=self.correction, correction_error=self.correction_error, T=self.T,
                            total_P=self.total_pressure)


if __name__ == "__main__":

    args = initialize()

    V = Viscosity(args.edr, args.gro, load=args.load)
    V.break_apart_trajectory(args.nsub)
    V.calculate()  # calculate viscosity
    if not args.load:
        V.save()  # save everything necessary to quickly reload data
    # V.plot_all()