#!/usr/bin/env python

import mdtraj as md
import argparse
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
from pymbar import timeseries
from scipy.optimize import curve_fit
import tqdm


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Trajectory file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-e', '--edr', default='wiggle.edr', type=str, help='Gromacs portable energy file')
    parser.add_argument('-T', '--temperature', default=300, type=float, help='Temperature of simulation')

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


def autocorrelation(x, largest_prime=500):
    """ FFT based autocorrelation function, which is faster than numpy.correlate. Efficiency is key in order to avoid
    headaches.
    :param x : numpy array of y values of an equispaced timeseries
    :param largest_prime : the largest prime factor of array length allowed. The smaller the faster. 1.6M points takes
    about 5 seconds with largest_prime=1000. Just be aware that you are losing data by truncating. But 5-6 data points
    isn't a big deal usually.
    """

    # Force length of array to be power of 2
    # powers = np.array([2 ** i for i in range(25)])
    # closest = np.argmin(np.abs(powers - (2*x.size - 1)))
    # nzeros = (powers[closest] // 2) - x.size
    # if nzeros < 0:
    #     closest += 1
    #     nzeros = (powers[closest] // 2) - x.size
    # v = np.concatenate((x, np.zeros([nzeros])))

    # Just don't allow a prime factor larger than 'largest_prime'. Truncate data until that condition is met
    l = 2*x.size - 1

    while largest_prime_factor(l) >= largest_prime or l % 2 == 0:
        l -= 1

    x = x[:(l + 1) // 2]
    length = x.size*2 - 1

    fftx = np.fft.fft(x, n=length)
    ret = np.fft.ifft(fftx * np.conjugate(fftx))
    ret = np.fft.fftshift(ret)

    auto = ret[length // 2:]  # used to return this, not sure if the following should happen
    n = ret[length // 2:].shape[0]
    average = [n - i for i in range(n)]

    return auto / average  # / np.amax(ret[length//2:])  # don't normalize


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


def running_integral(x):
    """
    :param x: numpy array of values at each time point
    :return: running integral of x.
    """

    integral = np.zeros_like(x)
    for i in range(x.shape[0]):
        integral[i] = np.sum(x[:i])

    return integral


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

    def __init__(self, edr, gro, max_lag=10):
        """
        :param edr: name of energy file (.edr)
        :param gro: name of coordinate file (.gro) (Assumes that system was run nve so volume is constant)
        :param max_lag : maximum time lag to for calculation of autocorrelation function (ps)
        """

        kb = 1.38064852 * 10 ** -23  # boltzmann constant
        t = md.load(gro)
        box = t.unitcell_vectors[0, :, :]
        self.volume = np.dot(box[2, :], np.cross(box[0, :], box[1, :]))  # volume of unit cell nm3
        self.volume *= 1 * 10 ** -27  # convert to m3
        self.max_lag = max_lag

        # get all off diagonal elements of the pressure tensor vs. time using gmx energy
        print('Extracting off-diagonal elements of the pressure tensor and Temperature from %s' % edr)
        ps = subprocess.Popen(('echo', 'Pres-XY', '\n', 'Pres-XZ', '\n', 'Pres-YX', '\n', 'Pres-YZ', '\n', 'Pres-ZX',
                               '\n', 'Pres-ZY', '\n', 'Temperature', '\n'), stdout=subprocess.PIPE)
        ps.wait()
        p = subprocess.Popen(('gmx', 'energy', '-f', '%s' % edr), stdin=ps.stdout,
                             stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        p.wait()

        with open('energy.xvg', 'r') as f:
            xvg = []
            for line in f:
                xvg.append(line)

        # find where the data starts by looking for the line that contains "Pres-ZY" since that was the last option
        # passed to gmx energy
        data_start = 0
        while xvg[data_start].count('Pres-ZY') == 0:
            data_start += 1
        data_start += 1

        P = np.zeros([len(xvg) - data_start, 6])
        T = np.zeros([len(xvg) - data_start])
        self.time = []
        for i in range(data_start, len(xvg)):
            data = xvg[i].split()
            self.time.append(float(data[0]))
            P[i - data_start, :] = [float(d) for d in data[2:]]
            T[i - data_start] = float(data[1])

        self.time = np.array(self.time)

        P *= 100000  # convert from bar to pascal
        # T_equil = timeseries.detectEquilibration(T)  # super slow for large number of data points
        self.Tavg = np.mean(T)
        self.prefactor = self.volume / (kb * self.Tavg)

        self.nintervals = 0
        while self.time[self.nintervals] < self.max_lag:
            self.nintervals += 1

        # autocorr = np.zeros_like(P)  # [=] Pa**2 * ps
        print('Calculating autocorrelation function for all 6 off-diagonals')
        # for i in tqdm.tqdm(range(6)):
        #     autocorr[:, i] = autocorrelation(P[:, i]).real

        self.autocorr = autocorrelation(P[:, 0]).real
        for i in tqdm.tqdm(range(1, 6)):
            self.autocorr += autocorrelation(P[:, i]).real  # calculate autocorrelation function w/FFT

        self.autocorr /= 6
        # self.autocorr = np.mean(autocorr, axis=1)  # average
        self.runningintegral = 1*10**-12 * running_integral(self.autocorr[:self.nintervals]) * self.prefactor

        bounds = ([-np.inf, 0, 0, 0], [np.inf, 1, np.inf, np.inf])  # bounds on parameters that will be fit

        solp, pcov = curve_fit(viscosity_fit, self.time[:self.nintervals], self.runningintegral[:self.nintervals],
                               bounds=bounds)
        self.fit_params = solp

        print('Infinite time viscosity: %.7f Pa-s' % viscosity_fit(np.inf, solp[0], solp[1], solp[2], solp[3]))

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
        #plt.ylabel('Running integral of autocorrelation function', fontsize=14)
        plt.legend()
        plt.tight_layout()
        if save:
            plt.savefig('running_integral.png')

        plt.show()


def autocorr(x):

    result = np.correlate(x, x, mode='full')

    return result[result.size//2:]


if __name__ == "__main__":

    args = initialize()

    V = Viscosity(args.edr, args.gro)
    V.plot_all()