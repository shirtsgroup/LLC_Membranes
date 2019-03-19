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
