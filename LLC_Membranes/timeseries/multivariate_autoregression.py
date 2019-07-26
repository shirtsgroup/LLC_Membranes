#!/usr/bin/env python

"""
Generate and/or fit multivariate autoregressive (MAR) processes
"""

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import timeseries


def generate_mar_process(phis, mu, cov, dim, ndraws):
    # TODO: move to generate_timeseries.py
    """ Create multivariate autoregressive timeseries. 'dim' dependent timeseries will be created with multivariate
    gaussian noise

    :param phis: autoregressive coefficients (order x dim x dim)
    :param mu: mean of trajectories in each dimension
    :param cov: covariance matrix describing multivariate normal distribution from which noise will be pulled
    :param dim: dimension, number of dependent trajectories to generate
    :param ndraws: number of points in each trajectory

    :type phis: numpy.ndarray
    :param mu: numpy.ndarray
    :type cov: numpy.ndarray
    :type dim: int
    :type ndraws: int
    """

    data = np.zeros([ndraws, dim])
    order = phis.shape[0]

    for d in range(order, ndraws):

        # calculate autoregressive terms
        data[d, :] = sum([phis[i, ...] @ data[d - (i + 1), :] for i in range(order)])

        # add gaussian noise
        data[d, :] += np.random.multivariate_normal(mu, cov)  # draw from multivariate normal distribution

    return data


def fit_mar_process(data, order):

    N, dim = data.shape

    Y = data[order:, :]  # (N - order)-by-d

    X = np.zeros([N - order, order, dim])
    for i in range(order, N):

        X[i - order, ...] = data[(i - order):i, :][:, ::-1]

    W = np.zeros([order, dim, dim])

    for i in range(dim):
        W[..., i] = np.linalg.inv(X[..., i].T @ X[..., i]) @ X[..., i].T @ Y

    return W


if __name__ == "__main__":

    m = 2  # order of autoregrssion (previous m terms are used to calculate value at current time)
    d = 2  # dimension. Number of dependent timeseries

    phis = np.zeros([m, d, d])  # autoregressive coefficients
    phis[0, ...] = np.array([[0.5, 0], [0, 0.4]])  # 1st lag
    phis[1, ...] = np.array([[0.3, 0], [0., -0.2]])  # 2nd lag

    mu = np.array([0, 0])  # mean of gaussian noise for each dimension
    cov = np.array([[1.5, 0.2], [0.2, 2]])  # covariance matrix for multivariate gaussian noise

    data = generate_mar_process(phis, mu, cov, d, 10000)

    vector_ar = timeseries.VectorAutoRegression(data, d)
