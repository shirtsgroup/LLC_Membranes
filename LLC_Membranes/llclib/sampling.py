#!/usr/bin/env python

from LLC_Membranes.llclib import fitting_functions
import numpy as np


def discrete_powerlaw_ccdf(val, xmin, alpha):
    """ Calculate the complementary cumulative distribution function of a discrete power
    law distribution such that P(x) = Pr(X >= x)

    :param x: value at which to evaluate CDF
    :param xmin: lower bound of power law PDF
    :param alpha: exponent of power law
    :param upper_limit: number of terms used to calculate Hurwitz zeta function

    :type x: float
    :type xmin: float
    :type alpha: float
    :type upper_limit: int

    :return CDF of power law PDF evaluated at x
    :rtype float
    """

    return fitting_functions.zeta(val, alpha) / fitting_functions.zeta(xmin, alpha)


def exact_discrete_power_law_sample(alpha, xmin, size=1):
    """ Exact method for generating random draws from a discrete power law distribution of
    form t**-alpha

    See Appendix D of https://epubs.siam.org/doi/abs/10.1137/070710111.

    This method is slow and the approximation might be more useful

    :param alpha: power law exponent
    :param xmin: lower limit of distribution.
    :param size: number of random draws to perform

    :type: alpha: float
    :type xmin: float
    :type size: int

    :return: array of random power law draws
    """

    r = np.random.uniform(0, 1, size=size)

    t = np.zeros([size])
    for i, val in enumerate(r):
        x2 = xmin
        while discrete_powerlaw_ccdf(x2, xmin, alpha) > (1 - val):
            x1 = x2
            x2 = 2 * x1
        t[i] = x2

    return t


def approximate_discrete_powerlaw(alpha, xmin, size=1):
    """ Approximate random draws from a discrete power law distribution of form t**-alpha.
    Much faster than discrete_powerlaw()

    See Appendix D of https://epubs.siam.org/doi/abs/10.1137/070710111

    :param alpha: power law exponent
    :param xmin: lower limit of distribution.
    :param size: number of random draws to perform

    :type: alpha: float
    :type xmin: float
    :type size: int

    :return: array of random power law draws
    """

    r = np.random.uniform(0, 1, size=size)

    t = np.round((xmin - 0.5) * (1 - r) ** (-1 / (alpha - 1)) + 0.5)

    return t


def random_power_law_dwell(alpha, ll=0.1, size=1, limit=None, discrete=False, exact=False):
    """ Randomly draw from a power law distribution of form t**-alpha
    See Appendix D of https://epubs.siam.org/doi/abs/10.1137/070710111

    :param alpha: anomalous exponent + 1
    :param ll: lower limit of distribution.
    :param size: number of random draws to perform
    :param limit: upper limit on dwell time length (As implemented, this is not mathematically sound, so it's only
    for demonstration purposes)
    :param discrete: pull from discrete power law distribution (uses approximation
    :param exact: use exact method for generating random draws from discrete power law distribution

    :type: alpha: float
    :type ll: float
    :type size: int

    :return: array of random power law draws
    """

    if limit is not None:

        r = 1 - np.exp((1 - alpha)*np.log(limit/ll))

    else:

        r = 1

    if discrete:

        if exact:

            return exact_discrete_power_law_sample(alpha, ll, size=size)

        else:

            return np.round((ll - 0.5)*(1 - np.random.uniform(0, r, size=size))**(-1/(alpha - 1)) + 0.5)

    else:

        return ll * (1 - np.random.uniform(0, r, size=size)) ** (-1 / (alpha - 1))


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