#!/usr/bin/env python

import numpy as np
from multiprocessing import Pool
import tqdm


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
    but uses FFTs. This method is faster than numpy.correlate. Efficiency is key in order to avoid headaches.

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
    autocorr_fxn /= np.arange(T.shape[0], 0, -1)[:, ...]

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

    size = len(x[0, :, axis].shape)  # number of axes in array where MSDs will be calculated

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
