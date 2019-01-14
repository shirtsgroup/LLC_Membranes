#! /usr/bin/env python

import numpy as np


def confidence_interval(data, confidence):
    """ Calculate confidence interval of data

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
