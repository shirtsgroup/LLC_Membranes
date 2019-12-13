#!/usr/bin/env python

import numpy as np
from multiprocessing import Pool
import tqdm
import warnings
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import stats, fitting_functions
from scipy.optimize import curve_fit


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


def autocovariance(t, largest_prime=500):

    return acf(t, largest_prime=largest_prime, autocov=True)


def acf_uneven(t, nboot=None, confidence=95.):
    """ Calculate autocorrelation function of uneven-length timeseries. Optionally generate bootstrapped statistics

    :param t: list of n time series. The length of each can be different.
    :param nboot: if not None, generate statistics using this many bootstrap trials
    :param confidence: confidence interval (expressed out of 100). Only necessary if nboot is not None

    :type t: list of lists
    :type nboot: None or int
    :type confidence: float

    :return: autocorrelation function and (optionally) bootstrapped statistics
    """

    len_ = np.array([len(x) for x in t])  # length of each trajectory

    max_len = max(len_)

    acf_ = np.zeros([len(t), max_len])

    keep = []  # list to hold indices of trajectories with a non-zero amount of hops
    for i in range(acf_.shape[0]):
        hops = t[i]
        if len(hops) > 2:  # correlation between two points is useless. Will always be +1, -1
            autocorrelation = acf(hops)
            acf_[i, :autocorrelation.size] = autocorrelation
            keep.append(i)

    acf_ = acf_[keep, :]
    len_ = len_[keep]

    if nboot is not None:

        boot = np.zeros([nboot, max_len])
        for b in range(nboot):
            sol = np.random.randint(acf_.shape[0], size=acf_.shape[0])
            for i in range(max(len_[sol])):
                ndx = sol[np.nonzero(acf_[sol, i])]
                if not list(ndx):
                    boot[b, i] = 0
                else:
                    #boot[b, i] = acf_[ndx, i].mean()
                    # give more weight to longer trajectories
                    boot[b, i] = np.average(acf_[ndx, i], weights=len_[ndx])
            # This is more pythonic, but multiple exceptions are raised so it doesn't really work as intended.
            # try:
            #     boot[b, i] = acf[ndx, i].mean()
            # except RuntimeWarning:  # happens if the solute with the max_hops dwell time is not included in 'sol'
            #     boot[b, i] = 0

        lower_confidence = (100 - confidence) / 2  # confidence intervals (percentiles)
        upper_confidence = 100 - lower_confidence

        errorbars = np.zeros([2, max_len])
        errorbars[0, :] = np.abs(np.percentile(boot, lower_confidence, axis=0) -
                                 boot.mean(axis=0))  # 2.5 percent of data below this value
        errorbars[1, :] = np.percentile(boot, upper_confidence, axis=0) - boot.mean(axis=0)

        return acf_, errorbars

    else:

        return acf_


def acf(t, largest_prime=500, autocov=False):

    """ Quickly calculated the autocorrelation function of a time series, t. This gives the same results as acf_slow()
    but uses FFTs. This method is faster than numpy.correlate.

    :param t: time series array (npoints, nseries)
    :param largest_prime : the largest prime factor of array length allowed. The smaller the faster. 1.6M points takes
    about 5 seconds with largest_prime=1000. Just be aware that you are losing data by truncating. But 5-6 data points
    isn't a big deal for large arrays.
    :param autocov: return autocovariance function insted (which is just the unnormalized autocorrelation)

    :type t: numpy.ndarray
    :type largest_prime: int
    :type autocov: bool

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

    if not autocov:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                autocorr_fxn /= np.var(T, axis=0)
            except FloatingPointError:
                print(autocorr_fxn)
                print(np.var(T, axis=0))
                exit()

    return autocorr_fxn  # normalized


def plot_autocorrelation(acfxn, errorbars=None, max_k=25, bootstrap=True, nboot=200, confidence=68.27, show=False,
                         fontsize=14, label=None, color='black', overlay=False):
    """ Plot autocorrelation function of increments

    :param acfxn: autocorrelation function of n trajectories (ntrajectories, npoints)
    :param errorbars: optional errorbars obtained by bootstrapping acf (2, npoints)
    :param max_k: maximum lag time to plot
    :param bootstrap: bootstrap data
    :param nboot: number of bootstrap trials
    :param confidence: confidence interval of shaded error region (percent)
    :param show: show the plot when done
    :param fontsize: size of font on axes and legend
    :param label: legend label
    :param color: line color
    :param overlay: If True, this will not create a new figure if one already exists.

    :type acfxn: numpy.ndarray
    :type max_k: int
    :type bootstrap: bool
    :type nboot: int
    :type confidence: float
    :type show: bool
    :type fontsize: int
    :type label: NoneType or str
    :type color: str
    :type overlay: bool
    """

    if not overlay:
        plt.figure()

    # calculate acf of each trajectory
    ntraj, n = acfxn.shape

    if bootstrap and errorbars is None:

        # This needs to be modified for uneven trajectory lengths
        boot = np.zeros([nboot, n])
        for i in range(nboot):
            ndx = np.random.randint(ntraj, size=ntraj)
            boot[i, :] = acfxn[ndx, :].mean(axis=0)
            #boot[i, :] = np.average(acfxn)

        errorbars = stats.confidence_interval(boot, confidence)

        avg = boot.mean(axis=0)
        plt.plot(np.arange(n), avg, lw=2, label=label, color=color)
        plt.fill_between(np.arange(n), avg + errorbars[0, :], avg - errorbars[1, :], alpha=0.25)

    else:

        plt.plot(np.arange(n), acfxn.mean(axis=0), lw=2, label=label, color=color)

        if errorbars is not None:
            plt.fill_between(np.arange(acfxn.shape[1]), errorbars[1, :] + acfxn.mean(axis=0), acfxn.mean(axis=0) -
                             errorbars[0, :], alpha=0.25)

    # formatting
    plt.xticks(np.arange(1, max_k)[::2])
    plt.xlim(-0.5, max_k)
    plt.ylim(-0.6, 1)
    plt.xlabel('Lag Time (steps)', fontsize=fontsize)
    plt.ylabel('Autocorrelation', fontsize=fontsize)
    plt.gcf().get_axes()[0].tick_params(labelsize=fontsize)
    plt.tight_layout()

    if show:
        plt.show()


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


def hurst(acf_, nboot=1, max_k=20, ll=1e-6):
    """ Estimate a distribution of Hurst parameters by fitting to an autocorrelation function. A distribution is
    obtained by bootstrapping.

    :param acf_: autocorrelation function for multiple trajectories (n x l) where n is the number of trajetories and
    l is the maximum length of the trajectories. The trajectories need not be the same length. This function can be
    obtained by timeseries.acf_uneven()
    :param nboot: number of bootstrap trials. This also determines the number of samples in the hurst distribution. If
    nboot=1 (default), you will get a single Hurst parameter based on a fit to the average acf.
    :param max_k: maximum time lag (in frames) used for fitting
    :param ll: lower limit on hurst parameter. If estimated H is below 0, then set the H equal to this lower limit. This
    should really only happen if you don't have a ton of data. And the dip in the acf should be pretty close to -0.5

    :type acf_: numpy.ndarray()
    :type nboot: int
    :type max_k: int
    :type ll: float

    :return: distribution of hurst parameters or a single hurst parameter if nboot=1
    :rtype: numpy.ndarray or float
    """

    weights = np.array([np.nonzero(acf_[i, :])[0].size for i in range(acf_.shape[0])])

    if nboot > 1:

        distribution = np.zeros([nboot])
        nsegments = acf_.shape[0]

        for b in range(nboot):

            traj = np.random.randint(0, nsegments, size=nsegments)
            #hboot = acf_[traj, 1].mean()  # first time lag autocovariance
            hboot = np.average(acf_[traj, 1], weights=weights[traj])
            H = np.log(2 * hboot + 2) / (2 * np.log(2))  # initial guess at H based on first dip in autocovariance

            acf_boot = np.zeros(max(weights[traj]))
            if H <= 0:  # a consequence of a small amount of data
                distribution[b] = ll
            else:
                distribution[b] = H
            # else:
            #     for i in range(max(weights[traj])):
            #         ndx = traj[np.nonzero(acf_[traj, i])]
            #         if not list(ndx):
            #             acf_boot[i] = 0
            #         else:
            #             # give more weight to longer trajectories
            #             acf_boot[i] = np.average(acf_[ndx, i], weights=weights[ndx])
            #
            #     #acf_boot = [acf_[traj, i][np.nonzero(acf_[traj, i])].mean() for i in range(max_k + 1)]
            #     if max_k > acf_boot.size:
            #         x = np.arange(acf_boot.size)
            #         y = acf_boot
            #     else:
            #         x = np.arange(max_k)
            #         y = acf_boot[:max_k]
            #
            #     h_opt = curve_fit(fitting_functions.hurst_autocovariance, x, y, p0=H)[0]
            #     distribution[b] = h_opt[0]

        return distribution

    else:

        # hboot = acf_[:, 1].mean()  # first time lag autocovariance
        hboot = np.average(acf_[:, 1], weights=weights)
        H = np.log(2 * hboot + 2) / (2 * np.log(2))  # initial guess at H based on first dip in autocovariance

        if H <= 0:
            h_opt = ll
        else:
            h_opt = curve_fit(fitting_functions.hurst_autocovariance, np.arange(max_k), acf_[:max_k], p0=H)[0]

        return h_opt[0]


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

        from statsmodels.tsa.api import VAR

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

    :type sequence: list or numpy.ndarray

    :return: list of indices where swithces between states occur
    :rtype: numpy.ndarray
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
