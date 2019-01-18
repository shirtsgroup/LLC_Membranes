#!/usr/bin/env python

import numpy as np
from scipy.stats import expon
from LLC_Membranes.analysis import Poly_fit
from sympy import mpmath


def log_power_law(x, a, alpha):
    """ Log power law dependence of MSD curve (for curve fitting)

    :param x: independent variable values
    :param a: scaling coefficient
    :param alpha: power (< 1: subdiffusive, 1: brownian, >1: superdiffusive)

    :return function evaluated at x values
    """

    return a + alpha * np.log(x)


def power_law(x, a, alpha):
    """ Power law dependence of MSD curve

    :param x: independent variable values
    :param a: scaling coefficient
    :param alpha: power (< 1: subdiffusive, 1: brownian, >1: superdiffusive)

    :return function evaluated at x values
    """

    return a * x ** alpha


def cdf_exp(left, right, scale):
    """ Calculate area under expnonential curve between two x locations """
    return expon.cdf(right, scale=scale) - expon.cdf(left, scale=scale)


def exponential_integrated(edges, A, B):
    """ Fit bins to exponential function based on integrated area under bins. The form of the exponential is: A*B*e^-Bx

    :param edges: edges of bins
    :param A: exponential function parameter. Controls scale of exponential function
    :param B: exponential function parameter. Controls rate of decay

    :type edges: np.ndarray
    :type A: float
    :type B: float

    :return integrated area between bin edges and under exponential curve as a function of x
    """

    bin_width = edges[1] - edges[0]  # bin width

    pdf = []
    for i in range(len(edges) - 1):
        pdf.append(cdf_exp(edges[i], edges[i + 1], 1/float(B)))

    pdf.append(1 - cdf_exp(0, edges[-1], 1/float(B)))

    return np.array(pdf) * (float(A) / bin_width)


def cdf_power_law(left, right, scale, alpha):
    """ Calculate area under power law curve between two x locations by integrating scale*e^-alpha from left to right.

    Scipy does not have tabulated data for power laws with alpha < 1, so I wrote my own for that case

    :param left: left bound of integral
    :param right: right bound of integral
    :param scale: parameter of power law function that controls its scale
    :param alpha: exponent of power law

    :type left: float
    :type right: float
    :type scale: float
    :type alpha: float

    :return: integrated area between left and right
    :rtype: float
    """

    exponent = 1 - alpha

    return (scale / exponent) * (right ** exponent - left ** exponent)


def powerlaw_integrated(edges, alpha, A):
    """ Fit bins to powerlaw function of form: At^-alpha

    :param edges: edges of bins to which power law will be fit
    :param alpha: exponent of power law
    :param A: parameter of power law function that controls its scale

    :type edges: np.ndarray
    :type alpha: float
    :type A: float

    :return: integrated area between bin edges and under power law curve as function of x
    """

    pdf = []
    for i in range(len(edges) - 1):
        pdf.append(cdf_power_law(edges[i], edges[i + 1], A, alpha))

    pdf.append(cdf_power_law(edges[-1], np.inf, A, alpha))

    return np.array(pdf) / (edges[1] - edges[0])


def fit_power_law(x, y, cut=1, interactive=True):
    """ Fit power law to MSD curves
    TODO: weighted fit (need to do error analysis first)

    :param y: y-axis values of MSD curve (x-axis values are values from self.time_uniform
    :param cut: fraction of trajectory to include in fit

    :type y: np.ndarray
    :type cut: float

    :return: Coefficient and exponent in power law of form [coefficient, power]
    """

    end = int(cut * len(x))  # fit up until a fraction, cut, of the trajectory

    nonzero = y.nonzero()

    # fit line to linear log plot
    A = Poly_fit.poly_fit(np.log(x[nonzero]), np.log(y[nonzero]), 1)[-1]

    return [np.exp(A[0]), A[1]]


def line(m, x, b):
    """ Return y values of a line given points, a slope and an intercept

    y = m*x + b

    :param m: slope
    :param x: points
    :param b: y-intercept

    :type m: float
    :type x: point or array of points
    :type b: float

    :return: y or array of y-values
    :rtype: np.ndarray()

    """

    return m * x + b


def zeta(x, alpha):
    """ Numerical estimate of Hurwitz zeta function. Used as discrete power law
    distribution normalization constant

    sum from n=0 to infinity of (n + x)^-\alpha

    :param x: point at which to evaluate function
    :param alpha: exponent of power law
    :param upper_limit: number of terms used to calculate Hurwitz zeta function (highest value of n above)

    :type x: float
    :type alpha: float
    :type upper_limit: int

    :return: evaluation of Hurwitz zeta function at x
    :rtype: float
    """

    if type(alpha) is np.ndarray:
        alpha = alpha[0]

    return float(mpmath.zeta(alpha, a=x))


def power_law_discrete_log_likelihood(alpha, x, xmin, minimize=False):
    """ Calculate log likelihood for alpha given a set of x values that might come from a
    power law distribution

    :param alpha: power law exponent. Calculates the log-likelihood of this value of alpha
    for the data
    :param x: array of values making up emperical distribution
    :param xmin: lower bound of power law distribution
    :param upper_limit: number of terms used to calculate zeta function

    :type alpha: float
    :type x: np.ndarray
    :type xmin: float
    :type upper_limit: int

    :return log-likelihood of input parameters
    :rtype float
    """

    n = x.size
    z = zeta(xmin, alpha)

    res = - n * np.log(z) - alpha * sum([np.log(i) for i in x])

    if minimize:
        return res * - 1
    else:
        return res


def gaussian_log_likelihood(parameters, data, maximize=False):
    """ Calculate log-likelihood given parameters and data

    :param parameters: a tuple of form (mean, sigma)
    :param data: data that might be gaussian
    :param maximize: if this is true, the opposite sign of the log-likelihood is returned so it can be used in a
    minimization function (as a way to calculate the maximum)

    :type parameters: tuple
    :type data: np.ndarray
    :type maximize: bool

    :return: log-likelihood
    """

    mean, sigma = parameters  # unpack parameters
    N = data.size

    L = -0.5 * N * np.log(2 * np.pi * sigma ** 2) - sum([((i - mean) ** 2) / (2 * sigma ** 2) for i in data])

    if maximize:
        return -L
    else:
        return L


def hurst_autocovariance(K, H):
    """ Return the analytical autocovariance of fractional gaussian noise for a given hurst exponent as a function
    of step number, k

    \gamma(k) = \dfrac{1}{2}[|k-1|^{2H} - 2|k|^{2H} + |k + 1|^{2H}]

    :param K: values of k at which to evaluate autocovariance (only integer values make sense)
    :param H: hurst exponent

    :return: analytical autocovariance function for fractional gaussian noise
    """

    return np.array([0.5 * (np.abs(k - 1)**(2*H) - 2 * np.abs(k)**(2*H) + np.abs(k + 1)**(2*H)) for k in K])
