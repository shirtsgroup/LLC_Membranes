#!/usr/bin/env python

from LLC_Membranes.llclib import timeseries
import numpy as np
import tqdm
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from statsmodels.tsa.stattools import adfuller


class Features(object):

    def __init__(self, xt):
        """

        :param xt: x-values of time series (nsteps, ntrajectories, ndimensions)

        """

        self.xt = xt
        self.ntraj = self.xt.shape[1]
        self.nframes = self.xt.shape[0]

        self.msd = None
        self.stationarity = None
        self.autocorrelation = None

    def compute_msd(self, ensemble=False, axis=(0, 1, 2), frac=0.4):
        """ Calculate means squared displacement of timeseries

        :param ensemble: calculate ensemble MSD instead of the default time-averaged MSD
        :param axis: which axes to include in MSD calculation (x = 0, y = 1, z = 2)

        :type ensemble: bool
        :type axis: int or tuple
        :type frac: time lag, expressed as a fraction of the simulation length, that will be reported as the MSD

        """

        self.msd = timeseries.msd(self.xt, axis, ensemble=ensemble)[int(frac*self.nframes)]

    def determine_stationarity(self, p=0.05):
        """ Use Augmented Dickey-Fuller test to determine if each time series is stationary. If p-value above p,
        timeseries is non-stationary. If p-value less than or equal to p, data is stationary.

        This function populates self.stationary with boolean values for each trajectory where True indicates a
        stationary process and False indicates a non-stationary process.

        :param p: p-value used as boundary for accepting / rejecting null hypothesis

        :type p: float

        :return:
        """

        self.stationarity = np.zeros([self.ntraj], dtype=bool)

        for t in tqdm.tqdm(range(self.xt.shape[1])):
            self.stationarity[t] = (adfuller(self.xt[:, t, 0])[1] < p)  # assumes axis=0 for now

    def compute_autocorrelation(self):
        """ Calculate autocorrelation between adjacent steps

        """

        self.autocorrelation = timeseries.acf(self.xt[1:, :, 0] - self.xt[:-1, :, 0])[1, :]

        # A noiseless CTRW trajectory during which no hops occur will have a variance in step size of zero. This will
        # result in an NaN during the autocorrelation function calculation. To remedy this, assume that there is some
        # white noise happening during the hop-less trajectory. The autocorrelation of white noise is 0.
        self.autocorrelation[np.isnan(self.autocorrelation)] = 0
