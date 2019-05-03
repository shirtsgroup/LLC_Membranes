#! /usr/bin/env python

import numpy as np
from scipy import stats


def confidence_interval(data, confidence):
    """ Calculate confidence interval of data.

    Plot these errorbars with plt.fill_between() as follows:
    plt.fill_between(x, mean + error[1, :], mean - error[0, :])


    :param data: array of data trajectories [n_trajectories, n_data_points]
    :param confidence: percent confidence

    :return: Upper and lower bounds to confidence intervals. Readily plotted with plt.errorbar
    """

    if type(data) is list:
        data = np.array(data)

    lower_confidence = (100 - confidence) / 2
    upper_confidence = 100 - lower_confidence

    mean_data = data.mean(axis=0)

    error = np.zeros([2, mean_data.size])  # [(lower,upper), number of data points
    error[0, :] = np.abs(np.percentile(data, lower_confidence, axis=0) - mean_data)  # percent of data below this value
    error[1, :] = np.percentile(data, upper_confidence, axis=0) - mean_data  # percent of data below this value

    return error


def outliers(data, alpha=0.01):
    """
    Check for outliers of viscosity calculation using Grubbs' test
    Steps:
    (1) Calculate critical t-statistic
    https://stackoverflow.com/questions/19339305/python-function-to-get-the-t-statistic?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    (2) Calculate critical G
    https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h1.htm

    :param data: data to search for outliers
    :param alpha : probability that point is falsely rejected (default=0.01)

    :type alpha: float
    :type data: numpy.ndarray

    :return: indices of outliers
    """

    outlier = True  # hypothesize that there is an outlier
    outliers = []
    x = np.copy(data)

    while outlier:
        n = len(x)  # number of samples
        t = stats.t.ppf(1 - ((alpha / 2) / (2 * n)), n - 2)
        gcrit = np.sqrt(t ** 2 / (n - 2 + t ** 2))
        gcrit *= (n - 1) / (n ** 0.5)
        G = np.abs((x - x.mean()) / x.std())
        potential_outlier = np.amax(G)
        ndx = np.argmax(G)
        if potential_outlier > gcrit:
            outliers.append(np.argmax(G))
            x = np.delete(x, ndx)
            #noutliers += 1
        else:
            outlier = False

    return outliers


class Cdf(object):  # taken from compare_disorder.py

    def __init__(self, data):
        """
        Generate an emperical cumulative distribution function from data
        :param data: x-values of data in no particular order
        """

        self.xs = np.array(sorted(data))
        self.N = float(len(self.xs))
        self.ys = np.arange(1, self.N + 1) / self.N

    def cdf(self, x):
        """
        Callable cumulative emperical distribution function
        :param x: array of x-values at which to evaluate cumulative emperical distribution function
        :return:
        """
        if type(x) is np.float64:
            x = np.array([x])

        ndx = [np.argmin(np.abs(self.xs - x[i])) for i in range(x.size)]

        return self.ys[ndx]

    def random_sample(self, n=1):
        """
        :param n: number of random samples to draw (default=1)
        :return: random samples
        """

        return np.random.choice(self.xs, size=n, replace=True)

    def update_cdf(self, obs):
        """ Add observation to the cdf (only works for single value right now)
        """

        self.xs = sorted(np.concatenate((self.xs, [obs])))
        self.N = float(len(self.xs))
        self.ys = np.arange(1, self.N + 1) / self.N
