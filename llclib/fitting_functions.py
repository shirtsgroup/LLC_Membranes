#!/usr/bin/env python

import numpy as np
from scipy.stats import expon
from LLC_Membranes.analysis import Poly_fit


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

    :return: Coefficient and exponent in power low of form [coefficient, power]
    """

    end = int(cut * len(x))  # fit up until a fraction, cut, of the trajectory

    # fit line to linear log plot
    A = Poly_fit.poly_fit(np.log(x[:end]), np.log(y[:end]), 1)[-1]

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
