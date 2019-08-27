#!/usr/bin/env python

import numpy as np
from multiprocessing import Pool
import tqdm
import warnings
import matplotlib.pyplot as plt
from statsmodels.tsa.api import VAR


def largest_prime_factor(n):
    i = 2
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
    return n


def acf_slow(d):
    """ Calculate the autocorrelation function of a time series. This speed of this method is O(n^2)

    :param d: numpy array of length n, with time series values {x1, x2 ... xn}
    :return: autocorrelation function
    """

    if type(d) is list:
        d = np.array(d)

    # Subtract mean
    d -= d.mean(axis=0)

    autocorr = np.zeros([len(d)])
    for l in range(d.shape[0]):  # cycle through lags
        N = d.shape[0] - l
        for n in range(N):
            autocorr[l] += d[n] * d[n + l]
        autocorr[l] /= N

    autocorr /= d.var()

    return autocorr


def acf(t, largest_prime=500):

    """ Quickly calculated the autocorrelation function of a time series, t. This gives the same results as acf_slow()
    but uses FFTs. This method is faster than numpy.correlate.

    :param t: time series : ndarray [npoints, nseries]
    :param largest_prime : the largest prime factor of array length allowed. The smaller the faster. 1.6M points takes
    about 5 seconds with largest_prime=1000. Just be aware that you are losing data by truncating. But 5-6 data points
    isn't a big deal for large arrays.

    """

    T = np.array(t)

    # Don't allow a prime factor larger than 'largest_prime'. Truncate data until that condition is met
    l = 2 * T.shape[0] - 1

    while largest_prime_factor(l) >= largest_prime or l % 2 == 0:
        l -= 1

    T = T[:(l + 1) // 2, ...]  # '...' allows for no second dimension if only a single time series is analysed
    length = T.shape[0] * 2 - 1

    T -= np.mean(T, axis=0)

    fftx = np.fft.fft(T, n=length, axis=0)
    ret = np.fft.ifft(fftx * np.conjugate(fftx), axis=0)
    ret = np.fft.fftshift(ret, axes=(0,))

    autocorr_fxn = ret[length // 2:].real

    if len(autocorr_fxn.shape) > 1:
        autocorr_fxn /= np.arange(T.shape[0], 0, -1)[:, None]
    else:
        autocorr_fxn /= np.arange(T.shape[0], 0, -1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        autocorr_fxn /= np.var(T, axis=0)

    return autocorr_fxn  # normalized


def autocov(joint_distribution, varied_length=False):

    """ Calculate the autocovariance function of the joint distribution of multiple realizations of a time series model

    See Pag 45 - 46 of Time Series Analysis (1st edition?) by James hamilton

    y_t : timeseries values at time t
    y_t-j : timeseries values at time t - j

    covariance_j = E(y_t - mu)(y_t-j - mu)

    In words: the covariance at lag j equals the expected value of y_t times y_t-j. They are not necessarily independent
    so you can't assume it equals E(y_t)*E(y_t-j)

    :param joint_distribution: n x m numpy array with n independent realizations of a time series consisting of m data
    points (observations) per realization.

    :returns autocovariance of joint distribution as function of lag j

    """

    observations = joint_distribution.shape[1]
    autocov = np.zeros([observations])
    counts = np.zeros([observations])

    for i in range(observations):

        nonzero = np.nonzero(joint_distribution[:, i])
        yt_expected_value = joint_distribution[nonzero, i] - joint_distribution[nonzero, i].mean()

        for j in range(observations):

            nonzero = np.nonzero(joint_distribution[:, j])
            ytlag_expected_value = joint_distribution[nonzero, j] - joint_distribution[nonzero, j].mean()

            autocov[abs(j - i)] += (yt_expected_value * ytlag_expected_value).mean()
            counts[abs(j - i)] += 1

    acov = autocov / counts

    # observations = joint_distribution.shape[1]
    # autocov = np.zeros([observations])
    # counts = np.zeros([observations])
    #
    # for i in range(observations):
    #
    #     yt_expected_value = joint_distribution[:, i] - joint_distribution[:, i].mean()
    #
    #     for j in range(observations):
    #
    #         ytlag_expected_value = joint_distribution[:, j] - joint_distribution[:, j].mean()
    #
    #         autocov[abs(j - i)] += (yt_expected_value * ytlag_expected_value).mean()
    #         counts[abs(j - i)] += 1
    #
    # acov = autocov / counts

    return acov / np.amax(acov)


def autocorrFFT(x):
    """ Function used for fast calculation of mean squared displacement

    :param x:
    :return:
    """

    N = len(x)
    F = np.fft.fft(x, n=2*N)  # 2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real   # now we have the autocorrelation in convention B
    n = N*np.ones(N) - np.arange(0, N)  # divide res(m) by (N-m)

    return res / n  # this is the autocorrelation in convention A


def msd_fft(args):
    """ Calculate msd using a fast fourier transform algorithm

    :param x: trajectory of particle positions, equispaced in time
    :param axis: axis along which to calculate msd ({x:0, y:1, z:2})

    :type x: np.ndarray
    :type axis: int

    :return: msd as a function of time
    """

    x, axis = args

    r = np.copy(x)
    r = r[:, axis]

    if len(r.shape) == 1:
        r = r[:, np.newaxis]

    N = len(r)
    D = np.square(r).sum(axis=1)
    D = np.append(D, 0)
    S2 = sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    Q = 2 * D.sum()
    S1 = np.zeros(N)
    for m in range(N):
      Q = Q - D[m - 1] - D[N - m]
      S1[m] = Q / (N - m)

    return S1 - 2 * S2


def msd_straightforward(x, axis):
    """
    Straightforward way to calculte msd. Gives same answer as msd()
    :param x: positions of centers of mass of all particles for each frame, numpy array [nframes, natoms, dim]
    :param ndx: list of indices to include in msd calculation (x = 0, y = 1, z = 2)

    :return: Average MSD and individual particle MSDs
    """
    n = x.shape[1]  # number of atoms
    N = x.shape[0]  # number of frames

    MSD = np.zeros([N])
    MSDs = np.zeros([N, n])

    for m in range(N):  # there nT different length intervals we can look at
        for k in range(N - m - 1):  # there are N - m - 1 independent intervals of length m
            MSDs[m, :] += np.linalg.norm(x[k + m, :, axis] - x[k, :, axis], axis=0)**2  # mean square displacement of all particles summed over all intervals of length m
        MSDs[m, :] /= (N - m)  # divide by the number of intervals to get an average msd over length m
        MSD[m] = np.mean(MSDs[m, :])

    return MSD, MSDs


def msd(x, axis, ensemble=False, nt=1):
    """ Calculate mean square displacement based on particle positions

    :param x: particle positions
    :param axis: axis along which you want MSD (0, 1, 2, [0, 1], [0, 2], [1, 2], [0, 1, 2])
    :param ensemble: if True, calculate the ensemble MSD instead of the time-averaged MSD

    :type x: ndarray (n_frames, n_particles, 3)
    :type axis: int or list of ints
    :type ensemble: bool

    :return: MSD of each particle
    """

    frames = x.shape[0]  # number of trajectory frames
    ntraj = x.shape[1]  # number of trajectories
    MSD = np.zeros([frames, ntraj], dtype=float)  # a set of MSDs per particle

    size = len(x[0, :, axis].shape)  # number of axes in array where MSDs will be calculate

    if ensemble:

        for n in range(ntraj):  # start at 1 since all row 0 will be all zeros
            MSD[:, n] = ensemble_msd(x[0, n, axis], x[:, n, axis], size)

    else:
        if nt > 1:
            with Pool(nt) as pool:
                for i, t in enumerate(pool.map(msd_fft, [(x[:, n, :], axis) for n in range(ntraj)])):
                    MSD[:, i] = t
        else:
            for n in tqdm.tqdm(range(ntraj)):
                MSD[:, n] = msd_fft((x[:, n, :], axis))

    return MSD


def ensemble_msd(x0, x, size):

    if size == 1:

        return (x - x0) ** 2

    else:

        return np.linalg.norm(x0 - x, axis=1) ** 2


def bootstrap_msd(msds, N, confidence=68):
    """ Estimate error at each point in the MSD curve using bootstrapping

    :param msds: mean squared discplacements to sample
    :param N: number of bootstrap trials
    :param confidence: percentile for error calculation

    :type msds: np.ndarray
    :type N: int
    :type confidence: float
    """

    nT, nparticles = msds.shape

    msd_average = msds.mean(axis=1)

    eMSDs = np.zeros([nT, N], dtype=float)  # create n bootstrapped trajectories

    print('Bootstrapping MSD curves...')
    for b in tqdm.tqdm(range(N)):
        indices = np.random.randint(0, nparticles, nparticles)  # randomly choose particles with replacement
        for n in range(nparticles):
            eMSDs[:, b] += msds[:, indices[n]]  # add the MSDs of a randomly selected particle
        eMSDs[:, b] /= nparticles  # average the MSDs

    lower_confidence = (100 - confidence) / 2
    upper_confidence = 100 - lower_confidence

    limits = np.zeros([2, nT], dtype=float)  # upper and lower bounds at each point along MSD curve
    # determine error bound for each tau (out of n MSD's, use that for the error bars)
    for t in range(nT):
        limits[0, t] = np.abs(np.percentile(eMSDs[t, :], lower_confidence) - msd_average[t])
        limits[1, t] = np.abs(np.percentile(eMSDs[t, :], upper_confidence) - msd_average[t])

    return limits


def step_autocorrelation(trajectories, axis=0):
    """ Calculate autocorrelation of step length and direction

    :param trajectories: array of position vs time (n_frames, n_particles, n_dimensions)
    :param axis: axis along which to calculate step lengths ({x:0, y:1, z:2})

    :type trajectories: numpy.ndarray
    :type axis: int or list
    """

    try:
        if len(axis) == 1:
            axis = axis[0]
    except TypeError:
        pass

    ntraj = trajectories.shape[1]  # number of particles with a trajectory

    # calculate acf of first trajectory in order to determine size of output array. timeseries.acf will truncate
    # the array slightly in order to make the FFT efficient
    ACF = acf(trajectories[1:, 0, axis] - trajectories[:-1, 0, axis])
    acfs = np.zeros([ntraj, ACF.size])
    acfs[0, :] = ACF

    keep = []
    for t in range(1, ntraj):
        steps = trajectories[1:, t, axis] - trajectories[:-1, t, axis]
        if not np.all(steps == 0):
            acfs[t, :] = acf(steps)
            keep.append(t)
        #acfs[t, :] = acf(trajectories[:ACF.size, t, axis])

    return acfs[keep, :]


def correlograms(zt):
    """ Plot correlograms of (z - zmean), (z - zmean)^2, (z - zmean)^3, (z - zmean)^4
    :param zt: timeseries of probability integral transforms

    :type zt: np.ndarray
    """

    # Initialize a 2x2 figure
    fig, ax = plt.subplots(2, 2, figsize=(9, 6), sharex=True, sharey=True)

    stop = 200  # index of point at which to stop (Diebold et al. stopped at 200)

    ax[0, 0].plot(acf(zt - zt.mean())[:stop], linewidth=2)  # (z - meanz)
    ax[0, 1].plot(acf((zt - zt.mean()) ** 2)[:stop], linewidth=2)  # (z - meanz)^2
    ax[1, 0].plot(acf((zt - zt.mean()) ** 3)[:stop], linewidth=2)  # (z - meanz)^2
    ax[1, 1].plot(acf((zt - zt.mean()) ** 4)[:stop], linewidth=2)  # (z - meanz)^2

    titles = ['$(z - \overline{z})$', '$(z - \overline{z})^2$', '$(z - \overline{z})^3$', '$(z - \overline{z})^4$']

    for i in range(4):
        ax[i // 2, i % 2].set_title(titles[i], fontsize=14)
        ax[i // 2, i % 2].tick_params(labelsize=14)

    ax[1, 0].set_xlabel('Lag (time steps)', fontsize=14)
    ax[1, 1].set_xlabel('Lag (time steps)', fontsize=14)
    ax[0, 0].set_ylabel('Correlation', fontsize=14)
    ax[1, 0].set_ylabel('Correlation', fontsize=14)


class VectorAutoRegression:

    def __init__(self, timeseries, r):
        r""" Fit a vector autogressive (VAR) process to data using statsmodels.tsa.vector_ar. The output object is
        just reduction and renaming of attributes produced after running the fit() method of the VAR class

        For more detailed docs, see: https://www.statsmodels.org/dev/vector_ar.html#module-statsmodels.tsa.vector_ar

        For a multidimensional time series, one could write a system of dependent autoregressive equations:

        .. math::

            Y_t = A_1*Y_{t-1} + ... + A_p*Y_{t-p} + u_t

        where

        .. math::

           Y_t = \begin{bmatrix} y_{1,t} \\ y_{2,t} \\ ... \\  y_{k,t} \end{bmatrix},
           Y_{t-1} = \begin{bmatrix} y_{1,t-1} \\ y_{2,t-1} \\ ... \\  y_{k,t-1} \end{bmatrix},
           ...

        The matrices :math:`A_i` are K x K matrices where K is the number of dimensions of the trajectory.
        :math:`A_1` contains the 1st time lag autoregressive coefficients. If

        .. math::

            A_1 = \begin{bmatrix} 0.5 & 0 \\ 0 & 0.4 \end{bmatrix}

        the associated system of equations for a VAR(1) process would be:

        .. math::

            y_{1,t} = 0.5y_{1,t-1} + u_{1,t}

            y_{2,t} = 0.4y_{2, t-1} + u_{2,t}

        Of course, adding cross-terms to A would create more complex dynamical behavior

        :math:`u_t` is a K-dimensional vector multivariate gaussian noise generated on the covariance matrix of the data

        :param timeseries: a T x K matrix where T is the number of observations and K is the number of variables/dimension
        :param r: autoregressive order. Number of past point on which current point depends

        :type timeseries: numpy.ndarray
        :type r: int
        """

        self.dim = timeseries.shape[1]  # number of dimensions

        # fit VAR model with statsmodels.tsa
        model = VAR(timeseries)
        results = model.fit(r)
        print(results.summary())

        # covariance matrix
        self.covariance = results.sigma_u_mle  # give same result as following commented out block
        # results summary stores the correlation matrix of residuals
        # https://blogs.sas.com/content/iml/2010/12/10/converting-between-correlation-and-covariance-matrices.html
        # corr = results.resid_corr  # residual correlation matrix
        # stds = results.resid.std(axis=0)  # standard deviation of data in each dimension
        # D = np.diag(stds)  # turn stds into a diagonal matrix
        # cov = D @ corr @ D  # convert correlation matrix to covariance matrix

        self.mu = results.params[0, :]
        self.mu_std = results.stderr[0, :]

        self.phi = results.coefs
        self.phi_std = np.zeros_like(self.phi)

        for i in range(self.dim):
            self.phi_std[:, i, :] = results.stderr[1:, i].reshape(r, self.dim)


def switch_points(sequence):
    """ Determine points in discrete state time series where switches between states occurs. NOTE: includes first and
    last point of time series

    :param sequence: series of discrete states

    :type sequence: list

    :return: list of indices where swithces between states occur
    :rtype: np.ndarray
    """

    # See https://stackoverflow.com/questions/36894822/how-do-i-identify-sequences-of-values-in-a-boolean-array
    switch_ndx = np.argwhere(np.diff(sequence)).squeeze().tolist()

    # add last frame as a switch point
    try:
        switch_ndx.append(len(sequence))
    except AttributeError:  # if there are no switches, it won't return a list
        switch_ndx = list([switch_ndx])
        switch_ndx.append(len(sequence))

    if switch_ndx[0] != 0:
        return np.array([0] + switch_ndx)  # also add first frame
    else:
        return np.array(switch_ndx)


def calculate_moving_average(series, n):
    """ Calculate moving average of a time series

    :param n: Number of previous points to average

    :type n: int
    """

    ret = np.cumsum(series, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]

    return ret[n - 1:] / n
