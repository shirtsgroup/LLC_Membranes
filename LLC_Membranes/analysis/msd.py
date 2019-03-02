#! /usr/bin/env python

import os
import sys
import argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from LLC_Membranes.analysis import Poly_fit, top
from LLC_Membranes.llclib import physical, topology, timeseries, fitting_functions, atom_props, file_rw
from scipy import stats
import tqdm
import sqlite3 as sql


def initialize():

    parser = argparse.ArgumentParser(description='Calculate mean squared displacement (MSD) and diffusion coefficient'
                                                 'for a specific residue or set of atoms in a trajectory')
    parser.add_argument('-t', '--trajectory', default='wiggle.trr', help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of .gro coordinate file')
    parser.add_argument('-r', '--residue', type=str, help='Name of residue whose diffusivity we want')
    parser.add_argument('-atoms', nargs='+', help='Name of atoms whose collective diffusivity is desired')
    parser.add_argument('-b', '--nboot', default=200, help='Number of bootstrap trials to be run')
    parser.add_argument('-f', '--frontfrac', default=0, type=float, help='Where to start fitting line on msd curve')
    parser.add_argument('-F', '--fracshow', default=.2, type=float, help='Percent of graph to show, also where to stop '
                        'fitting line during diffusivity calculation')
    parser.add_argument('-a', '--axis', default='z', type=str, help='Which axis to compute msd along')
    parser.add_argument('-nboot', default=200, type=int, help='Number of bootstrap trials for error estimation')
    parser.add_argument('-ensemble', '--ensemble', action="store_true", help='Calculate MSD as ensemble average')
    parser.add_argument('-compare', '--compare', action="store_true", help='Compare time-averaged and ensemble-averaged '
                                                                           'time series')
    parser.add_argument('-power_law', '--power_law', action="store_true", help='Fit MSD to a power law')

    parser.add_argument('-tails', '--tails', action="store_true", help='Only look at transport within tails')
    parser.add_argument('-pores', '--pores', action="store_true", help='Only look at transport within pores')
    parser.add_argument('-pr', '--pore_radius', default=1.5, type=float, help='Radius of pores. Anything greater than '
                        'this distance from the pore center will not be included in calculation')
    parser.add_argument('-acf', '--autocorrelation', action="store_true", help='Plot step autocorrelation function')
    parser.add_argument('-acov', '--autocovariance', action="store_true", help='Plot step autocovariance function')
    parser.add_argument('-nofit', '--nofit', action="store_true", help='Do not attempt to fit any curve to MSD')

    parser.add_argument('--update', action="store_true", help="update database with MD MSD values")
    parser.add_argument('-wt', '--wt_water', default=10, type=float, help='Weight percent water in system. Only need'
                                                                          'to adjust this if updating database.')
    parser.add_argument('-s', '--savename', default=False, type=str, help='If specified, save the MSD curve in a '
                                                                          'pickled .pl file of this name.')

    return parser


class Diffusivity(object):

    def __init__(self, traj, gro, axis, begin=0, startfit=0.01, endfit=1, residue=False, atoms=[], restrict=[]):
        """
        Calculate diffusivity from trajectory
        :param traj: unwrapped trajectory (i.e. gmx trjconv with -pbc nojump)
        :param gro: representative coordinate file
        :param axis: axis along which to compute MSD
        :param startfit: start linear fit to MSD startfit % into trajectory
        :param endfit: end linear fit to MSD endfit % into trajectory
        :param residue: if specified, the residue whose center of mass MSD will be measured
        :param atoms: if specified, group of atoms whose center of mass MSD will be measured
        :param restrict: restrict selection to certain indices. For example, if you want to calculate MSD of a certain
        residue, but only include a fraction of the total residues in the system.
        """

        # initialize path locations
        self.script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        self.top_location = "%s/../top/topologies" % self.script_location

        # initialize arguments for use in other functions
        self.gro = gro
        self.traj = traj
        self.residue = residue

        # initialize trajectory properties
        print('Loading trajectory...', end='', flush=True)
        self.t = md.load(self.traj, top=self.gro)[begin:]  # load trajectory
        print('Done!')
        self.nT = self.t.n_frames  # number of frames
        self.time = self.t.time / 1000  # time stamp on each frame, converted to nanoseconds

        # initialize fits to data, error analysis and plotting parameters
        self.startfit = int(startfit*self.nT)  # index at which to start fit
        if endfit == 1:
            self.endfit = self.nT - 1
        else:
            self.endfit = int(endfit*self.nT)  # index at which to end fit

        self.y_fit = []
        self.A = []  # fit parameters
        self.W = []  # weight matrix for curve fitting
        self.power_law_fit = None
        self.errorevery = int(np.ceil(self.nT / 100.0))  # plot only 100 bars total
        self.confidence_interval = 0
        self.yfit = None

        # autocorrelation / autocovariance
        self.acf = None
        self.acov = None

        # initialize results
        self.MSD = None
        self.limits = None
        self.MSD_average = 0
        self.slope_error = 0
        self.Davg = 0

        if residue:

            res = residue
            if res == 'SOL':  # mdtraj changes the name from SOL to HOH
                res = 'HOH'

            if restrict:
                selection = [a.index for a in self.t.topology.atoms if a.residue.name == res and a.index in restrict]
            else:
                selection = [a.index for a in self.t.topology.atoms if a.residue.name == res]

            topol = topology.Residue(res)
            atoms_per_residue = topol.natoms  # number atoms in a single residue
            matoms = [i for i in topol.mass.values()]  # mass of the atoms in residue
            self.mres = sum(matoms)

        elif atoms:

            selection = [a.index for a in self.t.topology.atoms if a.name in atoms]

            atoms_per_residue = len(atoms)
            matoms = np.array([atom_props.mass[x] for x in atoms])
            self.mres = np.sum(matoms)
        else:
            sys.exit('Error: No valid group of atoms or residues selected')

        self.map = topology.map_atoms(selection, nres_atoms=atoms_per_residue)
        pos = self.t.xyz[:, selection, :]

        self.axis = []
        if 'x' in axis:
            self.axis.append(0)
        if 'y' in axis:
            self.axis.append(1)
        if 'z' in axis:
            self.axis.append(2)

        print('Calculating center of mass of residues...', end='', flush=True)
        self.com = np.zeros([self.nT, pos.shape[1] // atoms_per_residue, 3])  # track the center of mass of each residue

        for f in range(self.nT):
            for i in range(self.com.shape[1]):
                w = (pos[f, i * atoms_per_residue:(i + 1) * atoms_per_residue, :].T * matoms).T  # weight each atom in the residue by its mass
                self.com[f, i, :] = np.sum(w, axis=0) / self.mres  # sum the coordinates and divide by the mass of the residue
        print('Done!')

        # plot 'random' z coordinate traces
        # np.random.seed(4)  # 4 gives a nice spread for ethanol
        # trajs = np.random.randint(0, self.com.shape[1], size=3)
        # print(trajs)
        # for i in trajs:
        #     plt.plot(self.time, self.com[:, i, 2], linewidth=2)
        #
        #     plt.ylabel('$z$-coordinate (nm)', fontsize=14)
        #     plt.xlabel('Time (ns)', fontsize=14)
        #     plt.gcf().get_axes()[0].tick_params(labelsize=14)
        # plt.tight_layout()
        #
        # plt.show()
        # exit()

        self.weights = True
        if self.com.shape[1] == 1:
            self.weights = False

        self.dt = self.time[-1] - self.time[-2]  # time step (assuming equispaced time points)

    def restrict_to_pore(self, r, dwell_fraction=0.95, tails=False, build_monomer='NAcarb11V.gro', spline=False,
                         buffer=0, npores=4):
        """ Restrict calculations to center of masses (COMs) that primarily stay in the pore OR tail region

        :param r: radius of pore. Anything greater than r from the pore center is considered the tail region
        :param dwell_fraction: Fraction of time spent in region of interest required in order to keep trajectory
        :param tails: if True, then restrict calculations to COMs primarily in the tail region
        :param build_monomer: monomer coordinate file of which liquid crystal membrane is mode
        :param spline: track pore centers with a 3D spline
        :param buffer: Do not count molecules below _buffer_ or above z-box-vector - _buffer_

        :type r: float
        :type tails: bool
        :type build_monomer: str
        :type spline: bool
        :type buffer: float
        """

        # find pore centers
        pore_defining_atoms = topology.LC(build_monomer).pore_defining_atoms
        pore_atoms = [a.index for a in self.t.topology.atoms if a.name in pore_defining_atoms]
        if spline:
            print('Creating pore splines')
            pore_centers = physical.trace_pores(self.t.xyz[:, pore_atoms, :], self.t.unitcell_vectors, 10)
        else:
            pore_centers = physical.avg_pore_loc(npores, self.t.xyz[:, pore_atoms, :], self.t.unitcell_vectors)

        inregion = physical.partition(self.com, pore_centers, r, buffer=buffer,
                                      unitcell=self.t.unitcell_vectors, npores=npores)

        if tails:
            inregion = ~inregion  # '~' flips True and False

        dwell = np.full((self.t.n_frames, self.com.shape[1]), False, dtype=bool)

        for t in range(self.t.n_frames):
            dwell[t, inregion[t]] = True

        fraction_dwelled = np.sum(dwell, axis=0) / self.t.n_frames  # fraction of total time spend in region of interest

        keep = np.where(fraction_dwelled >= dwell_fraction)[0]

        self.com = self.com[:, keep, :]

    def calculate(self, ensemble=False):

        print('Calculating MSD...', end='', flush=True)
        self.MSD = timeseries.msd(self.com, self.axis, ensemble=ensemble)
        self.MSD_average = np.mean(self.MSD, axis=1)
        print('Done!')

    def step_autocorrelation(self):
        """ Calculate autocorrelation of step length and direction
        """

        self.acf = timeseries.step_autocorrelation(self.com, axis=self.axis)

    def fit_linear(self):

        fit = 0
        while fit == 0:

            self.yfit, _, self.slope_error, _, A = Poly_fit.poly_fit(self.time[self.startfit:self.endfit],
                                                  self.MSD_average[self.startfit:self.endfit], 1, self.W)

            plt.plot(self.time[self.startfit:self.endfit], self.yfit, '--', color='black', label='Linear Fit')

            # plt.errorbar(self.time, self.MSD_average, yerr=[self.limits[0, :], self.limits[1, :]],
            #              errorevery=self.errorevery, label='MSD')
            plt.plot(self.time, self.MSD_average, label='MSD')

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
                self.startfit = int(self.startfit / (self.dt))  # convert time to index in t.time
                self.endfit = int(self.endfit / (self.dt))
            plt.clf()

    def fit_power_law(self, y):
        """ Fit power law to MSD curves
        :return: Coefficient and exponent in power low of form [coefficient, power]
        """

        A = Poly_fit.poly_fit(np.log(self.time[1:]), np.log(y[1:]), 1, self.W)[-1]

        return [np.exp(A[0]), A[1]]

    def bootstrap_power_law(self, N):

        self.power_law_fit = np.zeros([N, 2])

        print('Bootstrapping power law fit...')
        for b in tqdm.tqdm(range(N)):

            choices = np.random.randint(0, self.MSD.shape[1], size=self.MSD.shape[1])
            bootstrapped_msd = self.MSD[:, choices].mean(axis=1)
            self.power_law_fit[b, :] = self.fit_power_law(bootstrapped_msd)

            # plt.plot(self.time/1000, bootstrapped_msd)
            # plt.plot(self.time[1:]/1000, power_law(self.time[1:], self.power_law_fit[b, 0],
            #                                        self.power_law_fit[b, 1]))
            # plt.show()

    def bootstrap(self, N):
        """
        Estimate error at each point in the MSD curve using bootstrapping
        :param N: number of bootstrap trials
        """

        eMSDs = np.zeros([self.nT, N], dtype=float)  # create n bootstrapped trajectories

        print('Bootstrapping MSD curves...')
        for b in tqdm.tqdm(range(N)):
            indices = np.random.randint(0, self.com.shape[1], self.com.shape[1])  # randomly choose particles with replacement
            for n in range(self.com.shape[1]):
                eMSDs[:, b] += self.MSD[:, indices[n]]  # add the MSDs of a randomly selected particle
            eMSDs[:, b] /= self.com.shape[1]  # Divide every timestep by Nparticles -- average the MSDs

        confidence = 68  # percent confidence interval
        lower_confidence = (100 - confidence) / 2
        upper_confidence = 100 - lower_confidence

        self.limits = np.zeros([2, self.nT], dtype=float)  # upper and lower bounds at each point along MSD curve
        # determine error bound for each tau (out of n MSD's, use that for the error bars)
        for t in range(self.nT):
            self.limits[0, t] = np.abs(np.percentile(eMSDs[t, :], lower_confidence) - self.MSD_average[t])
            self.limits[1, t] = np.abs(np.percentile(eMSDs[t, :], upper_confidence) - self.MSD_average[t])

        npts = self.endfit - self.startfit
        if self.weights:
            self.W = np.zeros((npts, npts))
            for i in range(npts):
                self.W[i, i] = 1 / ((self.limits[0, i + self.startfit]) ** 2)
        else:
            self.W = 'none'

        slopes = np.zeros([N])
        for b in range(N):
            # fit line to each bootstrapped MSD
            A = Poly_fit.poly_fit(self.time[self.startfit:self.endfit], eMSDs[self.startfit:self.endfit, b],
                                  1, self.W)[-1]
            slopes[b] = A[1]

        slopes /= (2*100000*len(self.axis))  # nm^2 / ns = 1.0e-5 cm^2/s

        # calculate 95 % confidence interval
        self.confidence_interval = stats.t.interval(0.95, N - 1, loc=slopes.mean(), scale=stats.sem(slopes))
        self.Davg = slopes.mean()

    def plot(self, axis, fracshow=0.5, savedata=False, show=False):

        plt.figure()
        self.endfit = int(fracshow * self.time.size)

        plt.plot(self.time[:self.endfit], self.MSD_average[:self.endfit], label='MSD')
        plt.fill_between(self.time[:self.endfit], self.MSD_average[:self.endfit] + self.limits[0, :self.endfit],
                         self.MSD_average[:self.endfit] - self.limits[1, :self.endfit], alpha=0.7)

        last = self.MSD_average[(self.endfit - 1)]

        print('Final frame MSD: %.2f [%.2f, %.2f]' % (last, last - self.limits[1, self.endfit - 1],
                                                      last + self.limits[0, self.endfit - 1]))

        if savedata:

            np.savez_compressed('msd.npz', time=self.time[:self.endfit], msd=self.MSD_average[:self.endfit],
                                yerr=[self.limits[0, :self.endfit], self.limits[1, :self.endfit]])

        plt.ylabel('MSD ($nm^2$)', fontsize=14)
        plt.xlabel('time (ns)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()
        plt.savefig('Diffusivity_%s.pdf' % axis)

        if show:
            plt.show(block=True)

    def plot_power_law(self, bins=25):

        A = self.fit_power_law(self.MSD_average)

        fig, ax = plt.subplots(1, 2, figsize=(8, 4))

        ax[0].plot(self.time, self.MSD_average, linewidth=2)
        ax[0].plot(self.time, fitting_functions.power_law(self.time, A[0], A[1]), linewidth=2,
                   label=r'Power law fit to $Ae^{\alpha}$')
        ax[0].set_ylabel('Ensemble MSD ($nm^2/ns$)', fontsize=14)
        ax[0].set_xlabel('Time (ns)', fontsize=14)
        ax[0].xaxis.set_tick_params(labelsize=14)
        ax[0].yaxis.set_tick_params(labelsize=14)
        ax[0].legend()

        ax[1].hist(self.power_law_fit[:, 1], bins=bins, density=True)
        ax[1].set_title(r'Bootstrapped values of $\alpha$', fontsize=14)
        ax[1].set_ylabel('Probability', fontsize=14)
        ax[1].set_xlabel(r'$\alpha$', fontsize=14)
        ax[1].xaxis.set_tick_params(labelsize=14)

        plt.tick_params(labelsize=14)
        plt.tight_layout()
        # print(self.power_law_fit[:, 1].mean())
        # plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/transport/figures/msd_power_law.pdf')
        # plt.show()
        # exit()

    def plot_autocorrelation(self, show=True):
        """ Plot autocorrelation function

        :param show: show plot

        :type show: bool

        :return:
        """

        plt.figure()
        plt.plot(self.time[:self.acf.shape[1]], self.acf.mean(axis=0))
        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Autocovariance', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()

    def step_autocovariance(self):
        """ Calculate autocovariance of fractional gaussian noise in the trajectories (i.e. the step lengths)
        """

        self.acov = timeseries.autocov((self.com[1:, :, self.axis] - self.com[:-1, :, self.axis]).T[0, ...])

    def plot_autocovariance(self, show=True):
        """ Plot autocovariance function

        :param show: show plot

        :type show: bool
        """

        plt.figure()
        plt.plot(self.time[:-1], self.acov, color='blue', linewidth=3)
        plt.xlabel('Time Lag', fontsize=14)
        plt.ylabel('Autocovariance', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()

    def update_database(self, wt_water, file="../timeseries/msd.db", tablename="msd", ensemble=False):
        """ Update SQL database with information from this run

        :param wt_water: weight percent of water in system
        :param file: relative path (relative to directory where this script is stored) to database to be updated
        :param tablename: name of table being modified in database
        :param ensemble: True if ensemble MSD was calculated

        :type wt_water: float
        :type file: str
        :type tablename: str
        :type ensemble: bool
        """

        connection = sql.connect("%s/%s" % (self.script_location, file))
        crsr = connection.cursor()

        check_existence = "SELECT COUNT(1) FROM %s WHERE name = '%s' and wt_water = %.1f" % (tablename, self.residue,
                                                                                             wt_water)

        output = crsr.execute(check_existence).fetchall()

        msd = self.MSD_average[self.endfit - 1]  # usually I would do :self.endfit so subtracting 1 here makes sense
        msd_lower = msd - self.limits[1, self.endfit - 1]
        msd_upper = msd + self.limits[0, self.endfit - 1]

        if ensemble:
            data_labels = ['MD_MSD', 'MD_MSD_CI_lower', 'MD_MSD_CI_upper']
        else:
            data_labels = ['MD_TAMSD', 'MD_TAMSD_CI_lower', 'MD_TAMSD_CI_upper']

        if output[0][0] > 0:

            update_entry = "UPDATE %s SET %s = %.3f, %s = %.3f, %s = %.3f, sim_length = %.2f where name = '%s' and " \
                           "wt_water = %.1f" % \
                           (tablename, data_labels[0], msd, data_labels[1], msd_lower, data_labels[2], msd_upper,
                            self.time[-1], self.residue, wt_water)

            crsr.execute(update_entry)

        else:

            fill_new_entry = "INSERT INTO %s (name, %s, %s, %s, wt_water, sim_length) VALUES ('%s', %.3f, %.3f, %.3f, " \
                             "%.1f, %.2f)" % \
                             (tablename, data_labels[0], data_labels[1], data_labels[2], self.residue, msd, msd_lower,
                              msd_upper, wt_water, self.time[-1])

            crsr.execute(fill_new_entry)

        connection.commit()
        connection.close()


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.compare:

        show = False

        D_ensemble = Diffusivity(args.trajectory, args.gro, args.axis, residue=args.residue, atoms=args.atoms)

        if args.pores or args.tails:  # do this if solutes are restricted to tails or pores
            D_ensemble.restrict_to_pore(args.pore_radius, tails=args.tails)

        D_ensemble.calculate(ensemble=True)
        # D_ensemble.fit_linear()  # make sure diffusivity is being measured from linear region of the MSD curve
        D_ensemble.bootstrap(args.nboot)
        D_ensemble.plot(args.axis, fracshow=args.fracshow)
        print('D = %1.2e +/- %1.2e cm^2/s' % (D_ensemble.Davg, np.abs(D_ensemble.Davg -
                                                                     D_ensemble.confidence_interval[0])))
    else:
        show = True

    D = Diffusivity(args.trajectory, args.gro, args.axis, residue=args.residue, atoms=args.atoms)

    if args.pores or args.tails:  # do this if solutes are restricted to tails or pores
        D.restrict_to_pore(args.pore_radius, tails=args.tails)

    D.calculate(ensemble=args.ensemble)

    if args.power_law and not args.nofit:
        D.bootstrap_power_law(args.nboot)
        D.plot_power_law()
    else:
        if not args.compare and not args.nofit:
            D.fit_linear()  # make sure diffusivity is being measured from the linear region of the MSD curve

    if args.autocorrelation:
        D.step_autocorrelation()
        D.plot_autocorrelation()

    if args.autocovariance:
        D.step_autocovariance()
        D.plot_autocovariance()

    D.bootstrap(args.nboot)

    if args.savename:
        file_rw.save_object(D, '%s.pl' % args.savename)

    if args.ensemble:
        args.fracshow = 1  # same amount of statistics at each frame

    show = True
    D.plot(args.axis, fracshow=args.fracshow, show=show)

    if args.update:
        D.update_database(args.wt_water, ensemble=args.ensemble)

    if not args.power_law and not args.nofit:
        print('D = %1.2e +/- %1.2e cm^2/s' % (D.Davg, np.abs(D.Davg - D.confidence_interval[0])))

    if args.compare:

        plt.figure()
        end_frame = int(args.fracshow * D.time.size)

        plt.fill_between(D.time[:end_frame], D.MSD_average[:end_frame] + D.limits[0, :end_frame],
                         D.MSD_average[:end_frame] - D.limits[1, :end_frame], alpha=0.7, label='Time-averaged MSD')

        end_frame = int(args.fracshow * D.time.size)
        plt.fill_between(D_ensemble.time[:end_frame], D_ensemble.MSD_average[:end_frame] +
                         D_ensemble.limits[0, :end_frame], D_ensemble.MSD_average[:end_frame] -
                         D_ensemble.limits[1, :end_frame], alpha=0.7, label='Ensemble-averaged MSD')

        plt.ylabel('MSD ($nm^2$)', fontsize=14)
        plt.xlabel('time (ns)', fontsize=14)
        plt.legend(fontsize=14, loc=2)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.tight_layout()
        #plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/transport/figures/ethanol_msd_comparison.pdf')
        plt.show(block=True)
