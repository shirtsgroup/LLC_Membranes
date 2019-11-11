#!/usr/bin/env python

"""
Generate random variates of various distributions based on draws from a uniform distribution
See: https://www.cse.wustl.edu/~jain/books/ftp/ch5f_slides.pdf

Given the same random seed, MATLAB and numpy will generate the same draws from a uniform distribution. This isn't the
case for other distributions. Instead, it's possible to generate other distributions based on U(0, 1) draws. In turn,
one can generate the same pseudorandom numbers in MATLAB and numpy using these functions and MATLAB equivalent functions
"""

import numpy as np
from LLC_Membranes.llclib import fitting_functions
import math
np.set_printoptions(precision=4, suppress=True)


def randombeta(a, b, size=1):
    """ Generate random variates from beta distribution given shape parameters a, b
    """

    X = np.zeros([size])

    for i in range(size):
        X[i] = randombetavariate(a, b)

    return X


def randombetavariate(a, b):
    """ Generate single random variate from beta distribution given shape parameters a, b

    NOTE: only implemented for a < 0 and b < 0
    """

    if a < 1 and b < 1:
        x, y = [1, 1]  # make sure it enters while loop
        while x + y > 1:
            u1, u2 = np.random.uniform(size=2)
            x = u1 ** (1 / a)
            y = u2 ** (1 / b)

        return x / (x + y)

    # This doesn't work
    # elif isinstance(a, int) and isinstance(b, int):
    #
    #     u = sorted(np.random.uniform(size=(a + b + 1)))
    #
    #     return u[a - 1]


def randomexponential(a, size=1):
    """ Generate random variate from exponential distribution. Uses inverse CDF method

    pdf: f(x) = (1/a)e^(-x/a)

    :param a: scale paramter (a > 0)
    """

    u = np.random.uniform(size=size)

    return -a * np.log(u)


def randomgammaint(shape, scale=1.0, size=1):
    """ draw from random gamma distribution with integer shape parameter
    """

    U = np.random.uniform(size=(shape, size))
    X = -(1 / scale) * np.log(U).sum(axis=0)

    return X


def randomgamma(shape, scale=1.0, size=1, tol=1e-5):
    """ Generate random variates from a gamma distribution

    :param shape: Shape of gamma distribution (should be > 0)
    :param scale: Scale of the gamma distribution (should be > 0)
    :param size: number of variates to return

    :type shape: float or array_like of floats
    :type scale: float or array_like of floats
    :type size: int
    """

    if not isinstance(shape, (list, np.ndarray)):
        shape = [shape]

    rv = np.zeros([len(shape), size])

    for i, s in enumerate(shape):

        decimal, integer = math.modf(s)  # break up shape into integer and decimal

        if integer > tol:
            intgamma = randomgammaint(int(integer), scale=scale, size=size)
        else:
            intgamma = np.zeros([size])

        if decimal > tol:
            x = randombeta(decimal, 1 - decimal, size=size)
            y = randomexponential(1, size=size)
            decgamma = scale * x * y
        else:
            decgamma = np.zeros([size])

        rv[i, :] = intgamma + decgamma

    return rv


def randomnormal(mu, sigma, size=1):
    """ Generate random variate from the normal distribution using the Box Muller technique

    :return:
    """

    u = np.random.uniform(size=(2, size))

    x1 = mu + sigma*np.cos(2*np.pi*u[0, :]) * np.sqrt(-2*np.log(u[1, :]))
    #x2 = mu + sigma*np.sin(2*np.pi*u[0, :]) * np.sqrt(-2*np.log(u[1, :]))  # second independt Gaussian draw

    return x1


def randomwishart(a, d):

    sqrth = np.sqrt(0.5)
    norm = randomnormal(0, 1, d*d).reshape(d, d).T
    cholX = sqrth * np.triu(norm)
    i = np.arange(0, d)
    diag = [np.sqrt(randomgamma(g)) for g in a - i*0.5]
    for i in range(d):
        cholX[i, i] = diag[i]

    return cholX


def randombinomial(n, p):
    """ Generate random variates form a binomial(n, p) distribution

    :param n: number of trials
    :param p: probability of success
    """

    u = np.random.uniform(size=n)  # generate n U(0, 1) variates

    return sum(u < p)  # return number of variates less than p


def randomdirichlet(a):
    """ Python implementation of randdirichlet.m using randomgamma fucnction

    :param a: vector of weights (shape parameters to the gamma distribution)
    """

    x = randomgamma(a)
    x /= x.sum(axis=0)

    return x


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


def random_powerlaw(alpha, ll=0.1, size=1, limit=None, discrete=False, exact=False):
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


def random_exponential(lam, size=1, xmin=0):
    """ Randomly draw from an exponential distribution
    See Appendix D of https://epubs.siam.org/doi/abs/10.1137/070710111

    :param lam: rate of decay
    :param size: number of random draws to perform
    :param xmin: lowest value sampled

    :type lam: float
    :type size: int
    :type xmin: float

    :return: array of random draws
    """

    return xmin - np.log(1 - np.random.uniform(0, 1, size=size)) / lam


def random_powerlaw_cutoff(alpha, lam, xmin=1, size=1):
    """ Generate random variate from power law distribution with an exponential cut-off

    .. math::

        f(x) = x^{-\alpha}e^{-\lambda x}

    :param alpha: exponent describing decay of power law (higher alpha, faster decay)
    :param lam: exponent describing decay of exponential (higher lambda, faster decay)
    :param xmin: lowest value sample
    :param size: number of random variates to produce

    :type alpha: float
    :type lam: float
    :type xmin: float
    :type size: int

    :return: array of random draws
    """

    x = random_exponential(lam, size=size, xmin=xmin)
    p = (x / xmin) ** - alpha  # acceptance probability
    u = np.random.uniform(size=p.size)
    replace = np.where(u - p > 0)[0]
    while replace.size > 0:
        x[replace] = random_exponential(lam, size=replace.size, xmin=xmin)
        p[replace] = (x[replace] / xmin) ** - alpha
        u[replace] = np.random.uniform(size=replace.size)
        replace = np.where(u - p > 0)[0]  # rejected samples

    return x


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    alpha = 1.37
    lam = 0.00275
    rv = random_powerlaw_cutoff(alpha, lam, size=10000)

    fitting_functions.powerlaw_cutoff_mle(rv)

    x = np.linspace(1, 20, 1000)
    y = x**(-alpha)*np.exp(-lam*x)
    plt.plot(x, y, '--', lw=2, color='black')
    plt.plot(x, x**(-alpha), '--', lw=2, color='red')

    counts, bins = np.histogram(rv, bins=100, range=(0, 20))
    bin_width = bins[1] - bins[0]
    counts = [i / (counts.max() * y[0]) for i in counts]
    bin_centers = [i + bin_width/2 for i in bins[:-1]]
    plt.bar(bin_centers, counts, width=bin_width)

    plt.show()
    exit()

    np.random.seed(1)

    a = []
    for i in range(10000):
        a.append(randomdirichlet([10, 2])[0, 0])
    plt.hist(a, bins=50)
    plt.hist(np.random.beta(1, 2, size=10000), bins=50, alpha=0.5)
    plt.show()
    exit()
    plt.hist(np.random.binomial(100, 0.5, size=10000), bins=40, alpha=0.5, range=(30, 70))
    plt.show()

