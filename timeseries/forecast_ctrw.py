#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import ruptures
from scipy.optimize import curve_fit
from LLC_Membranes.llclib import topology, physical, file_rw, fitting_functions
from LLC_Membranes.analysis import Poly_fit, p2p
import tqdm
from scipy.stats import expon


def initialize():

    parser = argparse.ArgumentParser(description='Model a continuous time random walk')

    # MD trajectory control
    parser.add_argument('-t', '--trajectory', default='PR_nojump.xtc', help='Path to input file.')
    parser.add_argument('-g', '--gro', default='em.gro', help='Name of .gro coordinate file.')
    parser.add_argument('-b', '--begin', default=0, help='First trajectory frame used for analysis.')
    parser.add_argument('-e', '--end', default=-1, help='Last trajectory frame used for analysis.')
    parser.add_argument('-step', '--step', default=1, help='Include every "step" frames')
    parser.add_argument('-ma', '--moving_average', default=False, type=int, help='Calculate a moving average of the '
                        'center of mass coordinate')

    # loading saved objects
    parser.add_argument('-load', '--load', default=False, help='Specify name of pickled object to load')

    # restrict to pores parameters
    parser.add_argument('-restrict', '--restrict_to_pores', action='store_true', help='Only look at residues which'
                        'stay in the pore (based on last simulation frame)')
    parser.add_argument('--pore_defining_atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Atoms'
                        'used to define pore centers')
    parser.add_argument('--pore_defining_residue', default='HII', type=str, help='residue to which pore_defining atoms '
                                                                                 'belong')
    parser.add_argument('-radius', '--pore_radius', default=1.48, type=float, help='Defined pore radius (nm)')

    parser.add_argument('-r', '--residue', default='ETH', type=str, help='Name of residue whose diffusivity we want')
    parser.add_argument('-atoms', nargs='+', help='Name of atoms whose collective diffusivity is desired')
    parser.add_argument('-nbins', default=25, type=int, help='Number of bins to bin hop and dwell distributions into')

    # ctrw simulation
    parser.add_argument('-n', '--nhops', default=10000, type=int, help='Number of hops to perform')

    # bootstrapping
    parser.add_argument('-nboot', '--nboot', default=200, type=int, help='Number of bootstrap trials to be run')
    # parser.add_argument('-f', '--frontfrac', default=0.05, type=float, help='Where to start fitting line on msd curve')
    # parser.add_argument('-F', '--fracshow', default=.2, type=float, help='Percent of graph to show, also where to stop '
    #                     'fitting line during diffusivity calculation')
    # parser.add_argument('-a', '--axis', default='xyz', type=str, help='Which axis to compute msd along')
    # parser.add_argument('-nboot', default=200, type=int, help='Number of bootstrap trials for error estimation')
    # # | not implemented |
    # # v                 v
    # parser.add_argument('--restrict_to_pores', action="store_true", help='Only look at residue within pores of HII'
    #                                                                      'membrane')
    # parser.add_argument('-radius', default=1, type=float, help='Radius of pores. Anything greater than this distance'
    #                     'from the pore center will not be included in calculation')

    return parser


def autocorrFFT(x):

    N = len(x)
    F = np.fft.fft(x, n=2*N)  # 2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real   # now we have the autocorrelation in convention B
    n = N*np.ones(N)-np.arange(0, N)  # divide res(m) by (N-m)

    return res / n  # this is the autocorrelation in convention A


def msd_fft(x):

    N = len(x)
    D = np.square(x)
    D = np.append(D, 0)
    S2 = autocorrFFT(x)
    Q = 2*D.sum()
    S1 = np.zeros(N)
    for m in range(N):
      Q = Q-D[m-1]-D[N-m]
      S1[m] = Q / (N-m)

    return S1-2*S2


def gaussian(points, mean, sigma, amplitude):
    return (amplitude / np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(
        -(points - mean) ** 2 / (2 * sigma ** 2))


def wait_times(t, A, lam):

    # calculate ecdf

    return A*np.exp(-lam * t)


def wait_times_empirical(t, A):

    return A*np.exp(-lambda_empirical * t)


def wait_times_cdf(t, lam):

    return 1 - np.exp(-lam * t)


def random_dwell_time(lam):

        return -np.log(1 - np.random.uniform()) / lam


def cdf(left, right, scale):
    """ Calculate area under expnonential curve between two x locations """
    return expon.cdf(right, scale=scale) - expon.cdf(left, scale=scale)


def exponential(edges, A, B):

    bin_width = edges[1] - edges[0]  # bin width

    pdf = []
    for i in range(len(edges) - 1):
        pdf.append(cdf(edges[i], edges[i + 1], 1/float(B)))

    pdf.append(1 - cdf(0, edges[-1], 1/float(B)))

    return np.array(pdf) * (float(A) / bin_width)


def confidence_interval(data, confidence):
    """ Calculate confidence interval of data

    :param confidence: % confidence
    :return: error bars
    """

    lower_confidence = (100 - confidence) / 2
    upper_confidence = 100 - lower_confidence

    mean_data = data.mean(axis=0)

    error = np.zeros([2, mean_data.size])  # [(lower,upper), number of data points
    error[0, :] = np.abs(np.percentile(data, lower_confidence, axis=0) - mean_data)  # percent of data below this value
    error[1, :] = np.percentile(data, upper_confidence, axis=0) - mean_data  # percent of data below this value

    return error


class System(object):

    def __init__(self, traj, gro, res, start=0, end=-1, step=1, ma=False):
        """ Using an MD trajectory, calculate fit parameters for the dwell time distribution and hop length distribution

        :param traj: name of MD trajectory (GROMACS .xtc or .trr)
        :param gro: name of .gro (or .pdb) coordinate file associated with same topology as traj
        :param res: name of residue that is hopping around
        :param start: first frame of trajectory to include
        :param end: last frame of trajectory to include
        :param step: include every step frames in calculations
        :param ma: calculate a moving average on the center of mass coordinate traces

        :type traj: str
        :type gro: str
        :type res: str
        :type start: int
        :type end: int
        :type step: int
        """

        # load trjaectory
        print('Loading trajectory...', end='', flush=True)
        self.t = md.load(traj, top=gro)[start:end:step]
        print('Done!')

        self.time = self.t.time / 1000  # time in nanoseconds
        self.dt = self.time[1] - self.time[0]  # time step

        keep = [a.index for a in self.t.topology.atoms if a.residue.name == res]  # get indices of atoms making up residue of interest
        self.res_start = keep[0]  # index where residue starts

        self.residue = topology.Residue(res)
        self.nres = len(keep) // self.residue.natoms  # number of residues in system
        self.mass = [v for v in self.residue.mass.values()]  # mass of atoms making up residue
        self.pos = self.t.xyz[:, keep, :]  # positions of residue atoms

        print('Calculating centers of mass...', end='', flush=True)
        self.com = physical.center_of_mass(self.pos, self.mass)  # center of mass of residues
        print('Done!')

        if ma:
            self.calculate_moving_average(ma)

        # initialize for later
        self.pore_centers = None
        self.dwell_times = []  # distribution of dwell times
        self.tail_dwells = []
        self.hop_lengths = []
        self.lambda_distribution = []  # distribution of lambda for poisson process
        self.hop_sigma_distribution = []  # distribution of standard deviation of hop lengths

        self.partition = None  # array telling whether a solute is in the pores or tails. True if in pores, else False

    def plot_random_trajectories(self, n=1):
        """ Plot random center of mass trajectories and print indices of atoms in residue defining that COM

        :param n: Number of random trajectories to plot

        :type n: int
        """

        for i in np.random.randint(self.nres, size=n):
            start = self.res_start + len(self.mass)*i
            end = start + len(self.mass)
            print('Inidices %s through %s' % (start, end))
            plt.plot(self.com[:, i, 2], linewidth=2)
            plt.xlabel('Frame')
            plt.ylabel('Coordinate')
            plt.show()

    def calculate_moving_average(self, n):
        """ Calculate moving average of a time series

        :param n: Number of previous points to average

        :type n: int
        """

        nT = self.com.shape[0]
        ma = np.zeros([nT - n + 1, self.com.shape[1], 3])

        for i in range(self.com.shape[1]):
            for d in range(3):
                ret = np.cumsum(self.com[:, i, d], dtype=float)
                ret[n:] = ret[n:] - ret[:-n]
                ma[:, i, d] = ret[n - 1:] / n

        self.com = ma

    def calculate_solute_partition(self, r=1.5, spline=False, membrane_residue='HII', write_tcl=True):
        """ Determine whether each COM is in the tail or pore region as a function of time

        :param r: radial distance from pore center where pore region transitions to tail region
        :param spline: calculate pore centers as function of z using a spline in each pore
        :param membrane_residue: if using spline, give the name of the liquid crystal residue used to make the membrane

        :type r: float
        :type spline: bool
        :type membrane_residue: str
        """

        membrane = topology.LC('%s.gro' % membrane_residue)  # object w/ attributes of LC making up membrane
        keep = [a.index for a in self.t.topology.atoms if a.name in membrane.pore_defining_atoms and a.residue.name
                == membrane.name]
        self.pore_centers = physical.avg_pore_loc(4, self.t.xyz[:, keep, :], self.t.unitcell_vectors, buffer=0,
                                                  spline=spline, npts=10, progress=True, bins=False)

        self.partition = physical.partition(self.com, self.pore_centers, r, unitcell=self.t.unitcell_vectors,
                                            spline=spline)

        if write_tcl:
            self.write_tcl()

    def write_tcl(self, frame=-1, name='partition.tcl'):
        """ Write a .tcl script to view solute partition in VMD
        """

        pore_indices = [self.res_start + self.residue.natoms * i for i in np.where(self.partition[frame, :])[0]]
        tail_indices = [self.res_start + self.residue.natoms * i for i in np.where(self.partition[frame, :] == False)[0]]

        with open(name, 'w') as f:
            f.write('color Display Background white\n')
            f.write('mol addrep 0\n')
            f.write('mol modselect 0 0 index')
            for i in pore_indices:
                f.write(' %s' % i)
            f.write('\n')
            f.write('mol modcolor 0 0 ColorID 0\n')
            f.write('mol modstyle 0 0 CPK 2.0 0.3 12.0 12.0\n')
            f.write('mol addrep 0\n')
            f.write('mol modselect 1 0 index')
            for i in tail_indices:
                f.write(' %s' % i)
            f.write('\n')
            f.write('mol modstyle 1 0 CPK 2.0 0.3 12.0 12.0\n')
            f.write('mol modcolor 1 0 ColorID 1\n')

    def write_spline_coordinates(self, frame=-1):
        """ write spline coordinates to a .gro file
        """

        coordinates = np.reshape(self.pore_centers[frame, ...], (self.pore_centers.shape[1]*self.pore_centers.shape[2],
                                                                 3))
        file_rw.write_gro_pos(coordinates, 'spline.gro', name='K')

    def hops_and_dwells(self, penalty=1, nframes_dwell=10):
        """ Find breakpoints then assemble lists of dwell times and hop lengths. See documentation for Ruptures python
        package: http://ctruong.perso.math.cnrs.fr/ruptures-docs/build/html/index.html for more options that can be
        added

        :param penalty: penalty for cost function minimization
        :param nframes_dwell: number of frames that a solute needs to stay in the region of interest in order to be
        analyzed for hops.

        :type penalty: float
        :type nframes_dwell: int

        """

        # Handle initial and final dwell times since they don't have a true beginning or end. If the dwell times are
        # longer than some limit, then they will be included so that we don't miss out on dwells that last on the order
        # of the entire length of the trajectory
        initial_dwells = []  # length of time particle initially dwells before its first hop
        final_dwells = []  # length of time particle dwells up until end of trajectory

        for j in tqdm.tqdm(range(self.com.shape[1])):

            # Get indices where res jumps between regions
            # See https://stackoverflow.com/questions/36894822/how-do-i-identify-sequences-of-values-in-a-boolean-array
            switch_points = np.argwhere(np.diff(self.partition[:, j])).squeeze().tolist()
            try:
                switch_points.append(self.partition.shape[0])
            except AttributeError:
                switch_points = [switch_points]
                switch_points.append(self.partition.shape[0])

            # Analyze sub-trajectories where solute is in defined pore region
            begin = 0
            for i, end in enumerate(switch_points[::2]):  # only at switch points after which solute is in pores
                if (end - begin) > nframes_dwell:  # only include hops/dwells from trajectories long enough to analyze

                    # movement_3d = np.linalg.norm(self.com[begin:end, j, :], axis=1)  # magnitude of distance travelled
                    bp = ruptures.detection.Binseg(model='l2', min_size=1, jump=1).fit_predict(self.com[begin:end, j, :]
                                                                                               , pen=penalty)
                    # ruptures.display(self.com[begin:end, j, :], bp, figsize=(10, 6))
                    # plt.show()

                    initial_dwells.append(bp[0])
                    if len(bp) > 1:  # handle case where the only break point is at the end of the simulation
                        final_dwells.append(bp[-1] - bp[-2])

                    for k in range(len(bp) - 2):  # exclude first and last segments

                        #self.dwell_times.append(self.time[bp[k + 1]] - self.time[bp[k]])
                        self.dwell_times.append(bp[k + 1] - bp[k])  # discrete form

                        if k > 0:  # need two different z coordinates to get a hop length
                            self.hop_lengths.append(np.mean(self.com[bp[k]:bp[k + 1], j, 2]) -
                                                    np.mean(self.com[bp[k - 1]:bp[k], j, 2]))
                try:
                    begin = switch_points[2*i + 1]
                except IndexError:
                    pass

        self.tail_dwells = np.array(initial_dwells + final_dwells)

    def fit_distributions(self, nbins=25, plot=True, nboot=200):
        """ Fit curves to dwell time and hop length distributions

        :param nbins: number of bins in histograms
        :param plot: show the plots of the averaged distributions and bootstrapped parameter estimates
        :param nboot: number of bootstrap trials
        :param min_pl: minimum histogram value for power law (can't be zero)

        :return:
        """

        # Reconstruct bootstrapped distributions by randomly sampling from data
        # keep all bootstrapping data to generate error bars
        hist_dwell = np.zeros([nboot, nbins + 1])
        bins_dwell_centered = np.zeros([nboot, nbins + 1])
        p_dwell = [1, len(self.dwell_times) / (self.nres * self.time[-1])]
        hist_jump = np.zeros([nboot, nbins])
        bins_jump_centered = np.zeros([nboot, nbins])

        dwell_bin_width = 0
        hop_bin_width = 0
        dwell_mean = 0
        A_dwell = 0  # scaling factor for decaying expontential
        A_hop = 0

        #np.savez_compressed('dwells.npz', dwell_times=self.dwell_times)
        hist, edges = np.histogram(self.dwell_times, range=(1, 25), bins=25)
        normalized = hist / len(self.dwell_times)
        plt.bar(edges[:-1], normalized, 1, align='edge')

        plt.figure()
        plt.hist(self.dwell_times, bins=50)
        plt.show()
        exit()

        i = 0
        while i < nboot:

            # dwell times
            try:
                print('Bootstrapping ... Trial %s\r' % (i + 1), end='', flush=True)

                # Recreate dwell time distribution by randomly choosing from all dwell times with replacement
                dwell_times_boot = np.random.choice(self.dwell_times, size=len(self.dwell_times), replace=True)

                # Define edges of histogram
                edges = np.linspace(min(dwell_times_boot), max(dwell_times_boot), nbins + 1)

                # Add long dwell times randomly selected from beginning and end of trajectory
                tail_dwells = self.tail_dwells[np.where(self.tail_dwells > max(dwell_times_boot))[0]]
                dwell_times_boot = np.concatenate((np.random.choice(tail_dwells, size=len(tail_dwells), replace=True),
                                                  dwell_times_boot))

                # bin the bins
                bins = np.digitize(dwell_times_boot, edges)
                hist_dwell[i, :], bins_dwell = np.histogram(bins, nbins + 1)

                bin_width = edges[1] - edges[0]

                bins_dwell_centered[i, :] = np.array([edges[i] + bin_width / 2 for i in range(len(edges))])

                p_guess = [1.5, hist_dwell[i, np.argmin(np.abs(1 - edges))]]
                bounds = [[1, 0], [2, np.inf]]

                A, alpha = fitting_functions.fit_power_law(bins_dwell_centered[i, :], hist_dwell[i, :])

                # solp_dwell, cov_x = curve_fit(fitting_functions.powerlaw_integrated, edges, hist_dwell[i, :], p_guess,
                #                              bounds=bounds)

                # plot individual fit
                # alpha, A = solp_dwell
                # print(alpha, A)
                # alpha = 1.5
                # A = hist_dwell[i, 0] / (min(dwell_times_boot) + (bin_width/2))**-alpha
                #
                # pdf = fitting_functions.powerlaw_integrated(edges, alpha, A)
                #
                # plt.bar(edges, pdf, bin_width, align='edge')
                # print(hist_dwell[i, 0], A, min(dwell_times_boot), bin_width)

                plt.bar(edges, hist_dwell[i, :] / bin_width, bin_width, align='edge')
                plt.plot(bins_dwell_centered[i, :], A*bins_dwell_centered[i, :]**(-alpha), color='black')
                plt.show()
                exit()

                self.lambda_distribution.append(solp_dwell[1])

                # jump lengths
                hop_lengths_boot = np.random.choice(self.hop_lengths, size=len(self.hop_lengths), replace=True)
                hist_jump[i, :], bins_jump = np.histogram(hop_lengths_boot, bins=nbins)
                bins_jump_centered[i, :] = np.array([(bins_jump[i + 1] + bins_jump[i])/2 for i in
                                                     range(len(hist_jump[i, :]))])
                p_hops = [np.mean(hist_jump[i, :]), np.std(hist_jump[i, :]), 1]
                solp_hops, cov_x = curve_fit(gaussian, bins_jump_centered[i, :], hist_jump[i, :], p_hops)

                self.hop_sigma_distribution.append(np.abs(solp_hops[1]))  # sometime curve_fit finds negative answer

                # if np.abs(solp_hops[1]) > 2:
                #     plt.hist(hist_jump[i, :], bins=25)
                #     plt.show()

                # Things for plotting

                dwell_bin_width += bins_dwell_centered[i, 2] - bins_dwell_centered[i, 1]
                hop_bin_width += bins_jump_centered[i, 1] - bins_jump_centered[i, 0]
                A_dwell += solp_dwell[0]*solp_dwell[1]
                dwell_mean += solp_hops[0]
                A_hop += solp_hops[2]

                i += 1

            except RuntimeError:  # sometimes bootstrapping gives unfittable data
                continue

        dwell_bin_width /= nboot
        hop_bin_width /= nboot
        dwell_mean /= nboot
        A_dwell /= nboot
        A_hop /= nboot

        if plot:

            fig, ax = plt.subplots(2, 2, figsize=(12, 8))

            ax[0, 0].bar(bins_dwell_centered.mean(axis=0), hist_dwell.mean(axis=0),
                         width=dwell_bin_width, align='edge')
            ax[0, 0].set_xlabel('Dwell Time (ns)', fontsize=14)
            ax[0, 0].set_ylabel('Frequency', fontsize=14)
            ax[0, 0].plot(bins_dwell_centered.mean(axis=0), wait_times(bins_dwell_centered.mean(axis=0), A_dwell,
                    np.mean(self.lambda_distribution)), '--', color='black', label='$\lambda_{fit}$ = %.3f $\pm$ %.3f'
                    'ns$^{-1}$' % (np.mean(self.lambda_distribution), np.std(self.lambda_distribution)))
            ax[0, 0].legend(fontsize=12)

            ax[0, 1].hist(self.lambda_distribution, bins=nbins)
            ax[0, 1].set_xlabel('Bootstrapped $\lambda$ (ns$^{-1}$)', fontsize=14)
            ax[0, 1].set_ylabel('Frequency', fontsize=14)

            ax[1, 0].bar(bins_jump_centered.mean(axis=0), hist_jump.mean(axis=0), width=hop_bin_width)
            ax[1, 0].plot(bins_jump_centered.mean(axis=0), gaussian(bins_jump_centered.mean(axis=0), dwell_mean,
                    np.mean(self.hop_sigma_distribution), A_hop), '--', label='Gaussian fit\n $\sigma$=%.2f $\pm$ %.2f '
                    'nm' % (np.mean(self.hop_sigma_distribution), np.std(self.hop_sigma_distribution)), color='black')
            ax[1, 0].set_xlabel('Hop Length ($z$-direction, nm)', fontsize=14)
            ax[1, 0].set_ylabel('Frequency', fontsize=14)
            ax[1, 0].legend(fontsize=12)

            ax[1, 1].hist(self.hop_sigma_distribution, bins=nbins)
            ax[1, 1].set_xlabel('Bootstrapped $\sigma$ (nm)', fontsize=14)
            ax[1, 1].set_ylabel('Frequency', fontsize=14)

            plt.tight_layout()

            plt.show()


class CTRW(object):

    def __init__(self, nhops, dwell_lambdas, hop_sigmas, ntrials):

        self.dwell_lambdas = dwell_lambdas
        self.hop_sigmas = hop_sigmas
        self.nhops = nhops
        self.ntrials = ntrials

        self.trajectories = np.zeros([self.ntrials, self.nhops, 2])  # last dimension is [time, z_position]
        self.trajectory_hops = np.zeros([self.ntrials, 2*self.nhops - 1, 2])
        self.time_uniform = np.zeros([self.nhops])
        self.z_interpolated = np.zeros([self.ntrials, self.nhops])  # separate from time_uniform to save memory
        self.msd = np.zeros([self.ntrials, self.nhops])
        self.bootstrapped_msd = None
        self.yfit = None
        self.startfit = 0
        self.endfit = -1
        self.D = 0

    def generate_trajectories(self):

        print('Generating Trajectories...')
        for i in tqdm.tqdm(range(self.ntrials)):
            # constrain mean of hop-length distribution to be zero
            z_position = np.cumsum(np.random.normal(loc=0, scale=np.random.choice(self.hop_sigmas), size=self.nhops))
            self.trajectories[i, :, 1] = z_position - z_position[0]  # make initial z equal to 0

            time = np.zeros([self.nhops])
            lambda_trial = np.random.choice(self.dwell_lambdas)
            for j in range(1, self.nhops):  # make initial time equal to 0
                time[j] = random_dwell_time(lambda_trial)  # hop at random time intervals according to poisson process
            time = np.cumsum(time)
            self.trajectories[i, :, 0] = time

            self.trajectory_hops[i, 1::2, 0] = time[1:]
            self.trajectory_hops[i, 2::2, 0] = time[1:]

            self.trajectory_hops[i, ::2, 1] = self.trajectories[i, :, 1]
            self.trajectory_hops[i, 1:-1:2, 1] = self.trajectories[i, :-1, 1]
            self.trajectory_hops[i, -1, 1] = self.trajectories[i, -1, 1]

        print('Interpolating Trajectories...')
        # make uniform time intervals with the same interval for each simulated trajectory
        max_time = np.min(self.trajectories[:, -1, 0])
        self.time_uniform = np.linspace(0, max_time, self.nhops)
        for t in tqdm.tqdm(range(self.ntrials)):
            for i, x in enumerate(self.time_uniform):
                time_index = np.argmin(np.abs(x - self.trajectories[t, :, 0]))
                if x - self.trajectories[t, time_index, 0] < 0:
                    time_index -= 1
                self.z_interpolated[t, i] = self.trajectories[t, time_index, 1]

    def calculate_msd(self):
        """ Calculate the mean squared displacement for each CTRW trajectory

        """
        print('Computing mean squared displacement of each simulated trajectory')
        for t in tqdm.tqdm(range(self.ntrials)):
            self.msd[t, :] = msd_fft(self.z_interpolated[t, :])

    def plot_trajectory(self, n, show=False, save=True, savename='ctrw_trajectory.pdf'):
        """ Plot a CTRW trajectory

        :param n: Trajectory number
        """

        plt.figure()
        # plt.plot(self.trajectory_hops[n, :, 0] / 1000000, self.trajectory_hops[n, :, 1], linewidth=2)
        plt.plot(self.trajectories[n, :, 0] / 1000000, self.trajectories[n, :, 1], linewidth=2)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.xlabel('Time (ms)', fontsize=14)
        plt.ylabel('$z$-coordinate (nm)', fontsize=14)
        plt.tight_layout()

        if show:
            plt.show(block=True)
        if save:
            plt.savefig(savename)

    def fit_msd(self, start=0, end=-1):
        """Interactively fit a line to the mean squared displacement curve

        """

        self.startfit = start
        self.endfit = end
        dt = self.time_uniform[1] - self.time_uniform[0]

        plt.figure(2)
        fit = 0
        while fit == 0:

            self.yfit = Poly_fit.poly_fit(self.time_uniform[self.startfit:self.endfit],
                                                      self.msd.mean(axis=0)[self.startfit:self.endfit], 1)[0]

            plt.plot(self.time_uniform[self.startfit:self.endfit], self.yfit, '--', color='black', label='Linear Fit')
            plt.plot(self.time_uniform, self.msd.mean(axis=0), label='MSD')

            plt.ylabel('MSD ($nm^2$)', fontsize=14)
            plt.xlabel('time (ns)', fontsize=14)
            plt.gcf().get_axes()[0].tick_params(labelsize=14)
            plt.legend(loc=2)
            plt.tight_layout()
            plt.ion()
            plt.show()
            fit = int(input("Type '1' if the fit looks good: "))
            if fit != 1:
                print('Press enter to following prompts to leave as is')
                self.startfit = float(input("Time to start fit (ns): ") or self.startfit)
                self.endfit = float(input("Time to stop fit (ns): ") or self.endfit)
                self.startfit = int(self.startfit / dt)  # convert time to index in t.time
                self.endfit = int(self.endfit / dt)
                plt.clf()
            else:
                plt.close(2)

    def bootstrap_msd(self, nboot=200):

        # The average MSD is a collective property, so each bootstrap trial should be an average of self.ntrials
        # ranodomly reconstructed simulated trajectories
        self.bootstrapped_msd = np.zeros([nboot, self.msd.shape[1]])
        for i in range(nboot):
            indices = np.random.choice(self.msd.shape[0], size=self.msd.shape[0], replace=True)
            self.bootstrapped_msd[i, :] = self.msd[indices, :].mean(axis=0)

        slopes = []
        for i in range(nboot):
            A = Poly_fit.poly_fit(self.time_uniform[self.startfit:self.endfit],
                                  self.bootstrapped_msd[i, self.startfit:self.endfit], 1)[-1]
            slopes.append(A[1])

        self.D = [np.mean(slopes) / (2*1*10**9), np.std(slopes) / (2*1*10**9)]  # divide by dimension and converted to m^2/s

    def plot_msd(self, CI=95, nerrorbars=50, fracshow=0.5, save=True):
        """ Plot averaged mean squared displacement with error bars

        :param CI: confidence interval for error bars
        :param nerrorbars: show this many error bars
        :param fracshow: fraction of MSD to plot
        :param save: save the figure and the msd raw data

        :type CI: float
        :type nerrorbars: int
        :type fracshow: float between 0 and 1
        :type save: bool

        """

        plt.figure()
        error = confidence_interval(self.bootstrapped_msd, CI)
        plt.errorbar(self.time_uniform, self.msd.mean(axis=0), yerr=error, errorevery=self.msd.shape[1] // nerrorbars,
                     linewidth=2, elinewidth=2)
        plt.title('Diffusivity: %1.2e $\pm$ %1.2e m$^2$/s' % (self.D[0], self.D[1]))
        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Mean squared displacement (nm$^2$)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()
        plt.savefig('msd_ctrw.pdf')
        np.savez_compressed('msd.npz', msd=self.msd.mean(axis=0), error=error, time=self.time_uniform)
        plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.load:

        sys = file_rw.load_object(args.load)

    else:

        sys = System(args.trajectory, args.gro, args.residue, start=args.begin, end=args.end, step=args.step,
                     ma=args.moving_average)

        sys.calculate_solute_partition(spline=False, membrane_residue='HII')

        sys.hops_and_dwells()

        file_rw.save_object(sys, 'forecast_%s.pl' % args.residue)

    # sys.hops_and_dwells(penalty=1)
    # file_rw.save_object(sys, 'forecast_%s.pl' % args.residue)

    # histogram of time spent in pore region. Kind of interesting
    # plt.hist(np.sum(sys.partition, axis=0) / sys.partition.shape[0], bins=50)
    # plt.show()

    sys.fit_distributions(nbins=args.nbins, nboot=args.nboot)
    exit()

    random_walks = CTRW(args.nhops, sys.lambda_distribution, sys.hop_sigma_distribution, 200)
    random_walks.generate_trajectories()
    random_walks.calculate_msd()
    random_walks.fit_msd()
    random_walks.bootstrap_msd()
    random_walks.plot_msd()
    # for i in range(10):
    #     random_walks.plot_trajectory(i, show=False, save=True, savename='ctrw_trajectory_%d.pdf' % i)

    print('Diffusivity: %1.2e +/- %1.2e m$^2$/s' % (random_walks.D[0], random_walks.D[1]))
