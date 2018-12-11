#!/usr/bin/env python

import numpy as np


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
