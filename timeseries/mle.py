#!/usr/bin/env python

from scipy import stats
import numpy as np


def gaussian_log_likelihood(x, mu, sigma):
    """
    Calculate the log likelihood that a set of values, x, comes from a gaussian distribution given parameters for the
    mean and standard deviation
    :param x: int, float or numpy array of floats
    :param mu: mean
    :param sigma: standard deviation
    :return: log-likelihood
    """

    if isinstance(x, int) or isinstance(x, float):
        n = 1
    else:
        n = len(x)

    return -(n/2)*np.log(2*np.pi) - (n/2)*np.log(sigma**2) - (1 / (2*sigma**2))*np.sum(np.square(x - mu))


x = np.array([2, 3, 4, 5, 7, 8, 9, 10])

print(gaussian_log_likelihood(x, 5, 3))
print(gaussian_log_likelihood(x, 8, 3))
