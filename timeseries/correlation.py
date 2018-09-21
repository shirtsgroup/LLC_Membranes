#!/usr/bin/env python

import numpy as np


def largest_prime_factor(n):
    i = 2
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
    return n


def acf_slow(d):
    """ Calculate the autocorrelation function of a time series. This speed of this method is O(n^2)

    :param d: numpy array of length n, with time series values {x1, x2 ... xn}
    :return: autocorrelation function
    """

    # Subtract mean
    d -= d.mean(axis=0)

    autocorr = np.zeros([len(d)])
    for l in range(d.shape[0]):  # cycle through lags
        N = d.shape[0] - l
        for n in range(N):
            autocorr[l] += d[n] * d[n + l]
        autocorr[l] /= N

    autocorr /= d.var()

    return autocorr


def acf(t, largest_prime=500):

    """ Quickly calculated the autocorrelation function of a time series, t. This gives the same results as acf_slow()
    but uses FFTs. This method is faster than numpy.correlate. Efficiency is key in order to avoid headaches.

    :param t: time series : ndarray [npoints, nseries]
    :param largest_prime : the largest prime factor of array length allowed. The smaller the faster. 1.6M points takes
    about 5 seconds with largest_prime=1000. Just be aware that you are losing data by truncating. But 5-6 data points
    isn't a big deal for large arrays.

    """

    T = np.array(t)

    # Don't allow a prime factor larger than 'largest_prime'. Truncate data until that condition is met
    l = 2 * T.shape[0] - 1

    while largest_prime_factor(l) >= largest_prime or l % 2 == 0:
        l -= 1

    T = T[:(l + 1) // 2, ...]  # '...' allows for no second dimension if only a single time series is analysed
    length = T.shape[0] * 2 - 1

    T -= np.mean(T, axis=0)

    fftx = np.fft.fft(T, n=length, axis=0)
    ret = np.fft.ifft(fftx * np.conjugate(fftx), axis=0)
    ret = np.fft.fftshift(ret, axes=(0,))

    autocorr_fxn = ret[length // 2:].real
    autocorr_fxn /= np.arange(T.shape[0], 0, -1)[:, ...]
    autocorr_fxn /= np.var(T, axis=0)

    return autocorr_fxn  # normalized


def autocov(joint_distribution):

    """ Calculate the autocovariance function of the joint distribution of multiple realizations of a time series model

    See Pag 45 - 46 of Time Series Analysis (1st edition?) by James hamilton

    y_t : timeseries values at time t
    y_t-j : timeseries values at time t - j

    covariance_j = E(y_t - mu)(y_t-j - mu)

    In words: the covariance at lag j equals the expected value of y_t times y_t-j. They are not necessarily independent
    so you can't assume it equals E(y_t)*E(y_t-j)

    :param joint_distribution: n x m numpy array with n independent realizations of a time series consisting of m data
    points (observations) per realization.
    :returns autocovariance of joint distribution as function of lag j

    """

    observations = joint_distribution.shape[1]
    autocov = np.zeros([observations])

    for j in range(observations):

        yt_expected_value = joint_distribution[:, -1] - joint_distribution[:, -1].mean()
        ytlag_expected_value = joint_distribution[:, -j] - joint_distribution[:, -j].mean()

        autocov[j] = (yt_expected_value * ytlag_expected_value).mean()

    return autocov