#!/usr/bin/env python

""" Find the mean first passage time (MFPT) of a type of particle
"""

from LLC_Membranes.llclib import timeseries, file_rw, rand
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import tqdm
from multiprocessing import Pool
from scipy.sparse import csr_matrix as sparse_matrix
import warnings
from scipy.stats import norm
import levy
import fbm

with warnings.catch_warnings():
    warnings.filterwarnings('ignore', "Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient.")
    warnings.filterwarnings('ignore', "SparseEfficiencyWarning")


# TODO: should move these to llclib.rand
class Hops:

    def __init__(self, distribution, params):

        self.params = params

        if distribution.lower() in ['gaussian', 'normal', 'brownian']:
            self.name = 'gaussian'
        elif distribution.lower() in ['fbm']:
            self.name = 'fbm'
        elif distribution.lower() in ['flm']:
            self.name = 'flm'
        elif distribution.lower() in ['levy']:
            self.name = 'levy'

    def random_hop(self):

        if self.name == 'gaussian':
            return norm.rvs()
        elif self.name == 'levy':
            return levy.random(2, 0)

    def trajectory(self, length):

        if self.name == 'fbm':
            return self.params[0] * fbm.FBM(length, self.params[1], method="daviesharte").fbm() \
                   / ((1.0 / length) ** self.params[1])


class Dwell:

    def __init__(self, distribution, params):

        if distribution is None:
            self.name = None
        elif distribution.lower() in ['power']:
            self.name = 'power'
        elif distribution.lower() in ['power_cut']:
            self.name = 'power_cut'
        elif distribution.lower() in ['exponential']:
            self.name = 'exponential'

        self.params = params

    def random_dwell(self, size=1):

        if self.name is None:
            return np.ones(size)
        elif self.name == 'power':
            return rand.random_powerlaw(self.params[0], ll=self.params[1], size=size, limit=None,
                                        discrete=False, exact=False)
        elif self.name == 'power_cut':
            return rand.random_powerlaw_cutoff(self.params[0], self.params[1], xmin=self.params[2], size=size)
        elif self.name == 'exponential':
            return rand.random_exponential(self.params[0], size=size, xmin=self.params[1])


class Flux:

    def __init__(self, L, n, hop_dist='fbm', dwell_dist=None, dt=1., sigma=1, nbins=25, nt=1, save=True,
                 load=False):
        """ Generate concentration profiles

        :param L: length of 1D path (a pore for example) that particle follows (nm).
        :param n: number of particle trajectories to generate. It should be at least the number expected to be in a pore
        at steady state, but more is always better.
        :param hop_dist: hop length distribution
        :param dwell_dist:
        :param dt: timestep to use while contsructing trajectories
        :param sigma: width of hop distribution per unit time. The width of the hop distribution for the draws made each
        timestep is sqrt(dt) * sigma
        :param nbins: number of bins in concentration profile
        :param nt: number of threads to use when generating trajectories
        :param save: save trajectories to disk
        :param load: load trajectories saved to disk


        :type L: float
        :type n: int
        :type hop_dist: str
        :type dwell_dist: str or NoneType
        :type dt: float
        :type sigma: float
        :type nbins: int
        :type nt: int
        :type save: bool
        :type load: bool
        """

        # trajectory building parameters
        self.length = L
        self.ntraj = n
        self.sigma = sigma
        self.trajectories = []
        self.time = []
        self.dt = dt
        self.passage_times = []
        self.nt = nt

        H = 0.5
        params = (np.sqrt(self.dt) * sigma, H)
        if hop_dist.lower() in ['gaussian', 'normal', 'brownian']:
            self.hop_dist = Hops('gaussian', params)
        elif hop_dist.lower() in ['levy']:
            self.hop_dist = Hops('levy', params)
        elif hop_dist.lower() in ['fbm']:
            self.hop_dist = Hops('fbm', params)

        # TODO: these need to be passed
        params = [0.1, 1]
        if dwell_dist is None:
            self.dwell_dist = Dwell(None, params)
        elif dwell_dist.lower() in ['power']:
            self.dwell_dist = Dwell('power', params)
        elif dwell_dist.lower() in ['power_cut']:
            self.dwell_dist = Dwell('power_cut', params)
        elif dwell_dist.lower() in ['exponential']:
            self.dwell_dist = Dwell('exponential', params)

        print('Generating Trajectories...')
        self.generate_trajectories(save=save)

        # flux simulation arrays
        self.positions = None
        self.time_uniform = None
        self.flux_in = None
        self.flux_out = None
        self.pore_concentration = 0  # concentration at pore entrance
        self.dz = 0
        self.steps = 0
        self.nbins = nbins
        self.concentration = np.zeros([0, self.nbins])
        self.nparticles_in_pore = np.array([])

        # for animated plotting
        self.patches = None
        self.bins = None
        self.ax = None
        self.fig = None

    def _trajectory_realizations(self, ntraj, sigma):

        np.random.seed()  # need a different random seed for each thread

        trajectories = []
        times = []
        passage_times = []

        n = 0

        while len(passage_times) < ntraj:
            print('\r%d/%d trajectories' % (len(passage_times), ntraj), end='')

            if self.hop_dist.name == 'gaussian':

                walk = [0]

                # eventually move this to Hop class
                while 0 <= walk[-1] < self.length:

                    walk.append(walk[-1] + sigma * np.random.normal())

                if len(walk) > 2:

                    walk = np.array(walk)

                    trajectories.append(walk[1:].tolist())
                    time = np.cumsum(self.dwell_dist.random_dwell(
                        walk.size)) * self.dt  # TODO : not sure if multiplying by dt is right
                    times.append(time - time[0])

                    if walk[-1] >= self.length:

                        n += 1
                        passage_times.append(len(walk))

            elif self.hop_dist.name == 'fbm':

                walk = self.hop_dist.trajectory(2**16)

                if np.abs(walk).max() >= self.length:

                    # find when the cross-over occurs
                    # cross = timeseries.switch_points(np.abs(walk) >= self.length)
                    # walk = walk[:(cross[1] + 2)]

                    traj, time, ptimes = self._crossings(walk)

                    for i, t in enumerate(traj):

                        trajectories.append(t)
                        time.append(time[i])

                    passage_times += ptimes

        print('\n')

        return trajectories, times, passage_times

    def _recursive_trajectory_realizations(self, ntraj, sigma):
        """ A second way of generating trajectories. It's done by generating a single trajectory that eventually
        makes it the length of the pore. Then it is chopped up into segements based on whether the trajectory is
        positive or negative. Each of those segments is pulled back to the origin and the segmentation is repeated on
        the segment.

        :param ntraj:
        :param sigma:
        :return:
        """

        trajectories = []
        times = []

        for _ in tqdm.tqdm(range(ntraj), unit='trajectories'):

            walk = [0]
            while np.abs(walk[-1]) < self.length:

                walk.append(walk[-1] + sigma * np.random.normal())

            walk = np.array(walk)

            traj, time = self._crossings(walk)

            trajectories += traj
            times += time

        return trajectories, times

    def _crossings(self, x):

        # TODO: this needs to make sure trajectories stay within self.length
        time = np.ones(len(x))

        #crosses_boundary = timeseries.switch_points(np.logical_and(x > 0, x <= self.length))

        crosses_zero = timeseries.switch_points(x > 0)
        crosses_pore = timeseries.switch_points(np.abs(x) > self.length)

        crosses = np.unique(sorted(crosses_zero.tolist() + crosses_pore.tolist()))

        trajectories = []
        times = []
        passage_times = []

        for i, cross in enumerate(crosses[1:]):

            start = crosses[i] + 1
            if i == 0:
                start = 0

            end = cross + 1  # we want to include the last point

            # print(i, x.size, start, end)

            # traj = np.abs(x[start:end])

            traj = np.copy(x[start:end])

            try:
                if x[end] - x[start] > self.length:
                    passage_times.append(end - start)
            except IndexError:
                pass

            #
            # if traj.min() < -self.length:
            #     print(traj.min(), traj[-1], np.argmin(traj), len(traj))
            #     exit()
            #

            factor = np.abs(traj[0]) // (self.length / 2)
            if factor > 0:
                factor = ((factor - 1) // 2) + 1
                if traj[0] < 0:
                    traj += factor*self.length
                else:
                    traj -= factor*self.length

            # try:
            #     if np.abs(traj[0]) > self.length:
            #         if traj[0] < 0:
            #             traj += self.length
            #         else:
            #             traj -= self.length
            # except IndexError:
            #     pass

            # Recursion!
            #if timeseries.switch_points(np.logical_and(traj > 0, traj <= self.length)).size > 2:

            crosses_zero = timeseries.switch_points(x > 0)
            crosses_pore = timeseries.switch_points(np.abs(x) > self.length)

            crosses2 = np.unique(sorted(crosses_zero.tolist() + crosses_pore.tolist()))

            if crosses2.size > 2:

            #if timeseries.switch_points(traj > 0).size > 2:
                # plt.plot(traj)
                # plt.vlines(timeseries.switch_points(np.logical_and(traj > 0, traj <= self.length)), traj.min(), traj.max())
                # plt.show()
                # print('hello')
                # exit()

            # if timeseries.switch_points(traj > 0).size > 2:
            #     # print(traj)
            #     # print(timeseries.switch_points(np.logical_and(traj > 0, traj <= self.length)))
            #     # print('hello')
            #     # exit()
            # if timeseries.switch_points(x > 0).size > 2:
                traj_, time_, ptimes_ = self._crossings(traj)
                trajectories += traj_
                times += time_
                passage_times += ptimes_

            else:

                if traj.size > 1:

                    # if traj.size > 10:
                    #     print(x[start], x[end])
                    #     plt.plot(traj)
                    #     plt.plot(x[start:end])
                    #     plt.show()

                    trajectories += [np.abs(traj).tolist()]  # + sigma * np.random.normal())
                    times += [np.cumsum(time[start:end]) - 1]

        return trajectories, times, passage_times

    def _crossings_old(self, x):

        # TODO: this needs to make sure trajectories stay within self.length
        time = np.ones(len(x))

        # crosses_boundary = timeseries.switch_points(np.logical_and(x > 0, x <= self.length))
        crosses_boundary = timeseries.switch_points(x > 0)

        trajectories = []
        times = []

        for i, cross in enumerate(crosses_boundary[1:]):

            start = crosses_boundary[i] + 1
            if i == 0:
                start = 0

            end = cross + 1  # we want to include the last point

            # print(i, x.size, start, end)

            traj = np.abs(x[start:end])

            # traj = x[start:end]
            # if traj.min() < -self.length:
            #     print(traj.min(), traj[-1], np.argmin(traj), len(traj))
            #     exit()
            #
            # try:
            #     if np.abs(traj[end]) > self.length:
            #         if traj[end + 1] < 0:
            #             print('\n', x[end + 1])
            #             exit()
            #         traj -= self.length
            #
            # except IndexError:
            #     pass

            # Recursion!
            # if timeseries.switch_points(np.logical_and(traj > 0, traj <= self.length)).size > 2:
            if timeseries.switch_points(traj > 0).size > 2:
                # plt.plot(traj)
                # plt.vlines(timeseries.switch_points(np.logical_and(traj > 0, traj <= self.length)), traj.min(), traj.max())
                # plt.show()
                # print('hello')
                # exit()

                # if timeseries.switch_points(traj > 0).size > 2:
                #     # print(traj)
                #     # print(timeseries.switch_points(np.logical_and(traj > 0, traj <= self.length)))
                #     # print('hello')
                #     # exit()
                # if timeseries.switch_points(x > 0).size > 2:
                traj_, time_ = self._crossings(traj)
                trajectories += traj_
                times += time_

            else:

                if traj.size > 1:
                    trajectories += [np.abs(x[start:end])]  # + sigma * np.random.normal())
                    times += [np.cumsum(time[start:end]) - 1]

        return trajectories, times

    def generate_trajectories(self, save=True, savename='trajectories.npz'):
        """ Generate Brownian trajectories

        :param nt: number of threads
        """

        sigma = np.sqrt(self.dt) * self.sigma

        pool = Pool(processes=self.nt)

        trajectories_per_thread = self.ntraj // self.nt

        arguments = [(trajectories_per_thread, sigma) for _ in range(self.nt)]

        result = pool.starmap(self._trajectory_realizations, arguments)
        pool.close()

        passage_times = []
        for thread in range(self.nt):
            self.trajectories += result[thread][0]
            self.time += result[thread][1]
            passage_times += result[thread][2]

        self.passage_times = np.array(passage_times) * dt

        # print(self.passage_times.size)
        # print(np.unique(self.passage_times).size)
        # print(self.passage_times.mean())
        # plt.hist(self.passage_times, 25)
        # plt.show()
        # exit()
        # print(self.passage_times.mean())
        # plt.hist(self.passage_times, bins=50)
        # plt.show()
        # exit()

        if save:
            file_rw.save_object((self.trajectories, self.time), 'trajectories.pl')
            #np.savez_compressed(savename, trajectories=self.trajectories, time=self.time)

    def discretize_time(self):
        # TODO: increase bins in time_uniform (e.g. make time steps one tenth of dt)

        ntraj = len(self.trajectories)
        nperthread = ntraj // self.nt

        arguments = []
        for n in range(self.nt):
            if n < (self.nt - 1):
                arguments.append((n * nperthread, (n + 1) * nperthread))
            else:
                arguments.append((n * nperthread, (n + 1) * nperthread + (ntraj % self.nt)))

        pool = Pool(processes=self.nt)
        pool.starmap(self._discretize_group, arguments)
        pool.close()

    def _discretize_group(self, start, end):

        for i in tqdm.tqdm(range(start, end)):

            t = self.trajectories[i]

            last = np.argmin(np.abs(t[-1] - self.time_uniform))
            tu = self.time_uniform[:(last + 1)]

            # find the indices of the time point closest to the interpolated uniform time series
            time_index = np.argmin(np.abs(t[:, np.newaxis] - tu), axis=0)

            # make sure that jumps don't occur in interpolated trajectory until jumps occur in real trajectory
            time_index[np.where(tu - t[time_index] < 0)[0]] -= 1
            time_index[np.where(time_index < 0)] += 1  # Relax above restriction for first time point

            self.trajectories[i] = self.trajectories[i][time_index]
            self.time[i] = tu

            # self.visualize_discretization(self.trajectories[0], t)
            # self.visualize_discretization(self.trajectories[0][time_index], tu)
            # plt.show()

    @staticmethod
    def visualize_discretization(z, time):

        # for visualizing hops
        trajectory_hops = np.zeros([2 * len(time) - 1, 2])

        trajectory_hops[1::2, 0] = time[1:]
        trajectory_hops[2::2, 0] = time[1:]

        trajectory_hops[::2, 1] = z
        trajectory_hops[1:-1:2, 1] = z[:-1]
        trajectory_hops[-1, 1] = z[-1]
        plt.plot(trajectory_hops[:, 0], trajectory_hops[:, 1])

    def concentration_from_histogram(self):
        """ Generate a concentration profile by histogramming the atomic """

        data = []
        for t in self.trajectories:
            data += t#.tolist()

        file_rw.save_object(data, 'data.pl')

        print(self.nbins, self.length)
        plt.hist(data, self.nbins, range=(0, self.length))
        plt.show()

    def simulate_flux(self, pore_concentration=1, dz=0.1, steps=2000, measure_flux=False):
        """ Simulate a flux into the pore, starting new particles after a specified time
        """

        self.pore_concentration = pore_concentration

        longest = max([len(x) for x in self.trajectories])
        total_trajectories = len(self.trajectories)
        print('%d total trajectories generated' % total_trajectories)

        # discretize position versus time
        self.positions = sparse_matrix((self.ntraj, longest + steps * self.pore_concentration))

        self.dz = dz

        if self.dwell_dist.name is not None:
            print('Discretizing Trajectory...')
            n = 1
            self.time_uniform = np.arange(longest * n) * self.dt / n
            self.discretize_time()

        self.concentration_from_histogram()
        exit()

        self.steps = steps

        if measure_flux:
            self.flux_out = np.zeros(self.steps)

        traj_no = 0
        start = 0
        naddtot = []
        inlet = np.zeros(steps)
        dwell_times = []

        print('Simulating flux by holding the interface concentration constant at %d particles' % self.pore_concentration)
        for step in tqdm.tqdm(range(steps), unit=' Time Steps'):

            nadd = 0
            if step > 0:
                #print(inlet[:step], inlet[:step].mean())
                if inlet[:step].mean() < pore_concentration:
                    nadd = 1  # This will need to be changed to handle higher inlet concentration
            else:
                nadd = 1

            # figure out how many solutes to add in order to keep concentration constant
            #if self._get_inlet_concentration(start) == 0 and self._get_inlet_concentration(start + 1) == 0:
            #    nadd = 1
            #else:
            #    nadd = 0

            
            #c = self._get_inlet_concentration(start)
            #nadd = self.pore_concentration - c
            #inlet.append(c)

            #if nadd < 0:  # pore concentration too high. Don't add any
            #    nadd = 0

            naddtot.append(nadd)

            if traj_no + nadd > self.ntraj:
                print("Cleaning position matrix ...")

                previous_step = self.nparticles_in_pore.size
                # Record the total number of particles in the pore before modifying position matrix
                self._update_nparticles_in_pore(step - previous_step)

                # Record concentration profile in pore before modifying position matrix
                self._update_concentration_profile(step - previous_step)

                # Remove trajetories that are already finished
                self.positions, traj_no = self._clean_position_matrix(step - previous_step)
                start = 0

            for particle in range(traj_no, traj_no + nadd):

                random_trajectory_index = np.random.choice(total_trajectories)  # randomly choose a trajectory

                traj = self.trajectories[random_trajectory_index]
                if traj[-1] >= self.length:
                    dwell_times.append(len(traj))
                # print(traj)
                # exit()
                time = self.time[random_trajectory_index]
                last = np.argmin(np.abs(time[-1] - self.time_uniform))
                # tu = self.time_uniform[:(last + 1)]

                # find the indices of the time point closest to the interpolated uniform time series
                # time_index = np.argmin(np.abs(time[:, np.newaxis] - tu), axis=0)
                # # make sure that jumps don't occur in interpolated trajectory until jumps occur in real trajectory
                # time_index[np.where(tu - time[time_index] < 0)[0]] -= 1
                # time_index[np.where(time_index < 0)] += 1  # Relax above restriction for first time point
                #
                # interpolated = traj[time_index]  # TODO: this has issues for trajectories starting at zero

                #self.positions[particle, start:interpolated.size + start] = interpolated
                self.positions[particle, start:(traj.size + start)] = traj

                # print(self._get_inlet_concentration(start))

            if measure_flux:

                self.flux_out[step] = self._count_flux_out(start)

            inlet[step] = self._get_inlet_concentration(start)
            #print(inlet[step], nadd)
            #data = self.positions[:, start].data
            #print(data)
            #print(np.where(np.logical_and(data >=0, data < self.dz))[0].size)
            #if step == 10:
            #    exit()

            traj_no += nadd
            start += 1

        print(naddtot[equil:])
        print(np.mean(naddtot[equil:]))

        print(sum(naddtot) / (steps * self.dt))
        plt.plot(naddtot)
        #plt.plot(self.flux_out)
        plt.figure()
        plt.plot(inlet)
        plt.title('Inlet Concentration')
        file_rw.save_object(inlet, 'inlet.pl')
        file_rw.save_object(dwell_times, 'dwell_times.pl')

        plt.figure()
        plt.hist(dwell_times, bins=25)

        previous_step = self.nparticles_in_pore.size
        self._update_nparticles_in_pore(steps - previous_step)
        self._update_concentration_profile(steps - previous_step)

        #self.animate_concentration_profile()

    def simulate_flux_constant_particles(self, n=50, steps=2000, measure_flux=False):
        """ Simulate a flux into the pore, starting new particles after a specified time
        """

        longest = max([len(x) for x in self.trajectories])
        total_trajectories = len(self.trajectories)

        # discretize position versus time
        self.positions = sparse_matrix((self.ntraj, longest + steps))

        self.time_uniform = np.arange(self.positions.shape[1]) * self.dt
        self.steps = steps

        if measure_flux:
            self.flux_out = np.zeros(self.steps)

        traj_no = 0
        start = 0
        naddtot = []

        dwell_times = []
        print('Simulating flux by holding the nubmer of particles constant at %d' % n)
        for step in tqdm.tqdm(range(steps), unit=' Time Steps'):

            if step > 0:
                nadd = n - self._get_nparticles_in_pore(start)
                if nadd < 0:
                    print(nadd)
                    exit()
            else:
                nadd = n

            naddtot.append(nadd)

            if traj_no + nadd > self.ntraj:
                print("Cleaning position matrix ...")

                previous_step = self.nparticles_in_pore.size
                # Record the total number of particles in the pore before modifying position matrix
                self._update_nparticles_in_pore(step - previous_step)

                # Record concentration profile in pore before modifying position matrix
                self._update_concentration_profile(step - previous_step)

                # Remove trajetories that are already finished
                self.positions, traj_no = self._clean_position_matrix(step - previous_step)
                start = 0

            for particle in range(traj_no, traj_no + nadd):

                random_trajectory_index = np.random.choice(total_trajectories)  # randomly choose a trajectory

                traj = self.trajectories[random_trajectory_index]
                if traj[-1] >= self.length:
                    dwell_times.append(len(traj))

                time = self.time[random_trajectory_index]
                last = np.argmin(np.abs(time[-1] - self.time_uniform))
                # tu = self.time_uniform[:(last + 1)]

                # find the indices of the time point closest to the interpolated uniform time series
                # time_index = np.argmin(np.abs(time[:, np.newaxis] - tu), axis=0)
                # # make sure that jumps don't occur in interpolated trajectory until jumps occur in real trajectory
                # time_index[np.where(tu - time[time_index] < 0)[0]] -= 1
                # time_index[np.where(time_index < 0)] += 1  # Relax above restriction for first time point
                #
                # interpolated = traj[time_index]  # TODO: this has issues for trajectories starting at zero

                # self.positions[particle, start:interpolated.size + start] = interpolated
                self.positions[particle, start:(traj.size + start)] = traj

                # print(self._get_inlet_concentration(start))

            if measure_flux:
                self.flux_out[step] = self._count_flux_out(start)

            traj_no += nadd
            start += 1

        file_rw.save_object(self.flux_out, 'flux.pl')

        plt.plot(naddtot)

        file_rw.save_object(dwell_times, 'dwell_times.pl')

        plt.figure()
        plt.hist(dwell_times, bins=25)

        previous_step = self.nparticles_in_pore.size
        self._update_nparticles_in_pore(steps - previous_step)
        self._update_concentration_profile(steps - previous_step)

    def _get_nparticles_in_pore(self, step):

        data = self.positions[:, step].data
        return np.where(np.logical_and(data >= 0, data < self.length))[0].size

    def _get_inlet_concentration(self, step):
        """ Number of particles in inlet slab

        :param step: step number

        :type step: int

        :return: number of particles in inlet slab
        """

        # this gets progressively slower as sparse matrix grows
        # return np.where(self.positions.tocsr()[:, step].data < self.dz)[0].size

        # for lil_matrix
        # return np.where(np.array(self.positions[:, step].T.data[0]) < self.dz)[0].size

        # when self.positions is already a csr matrix
        #return np.where(self.positions[:, step].data < self.dz)[0].size
        data = self.positions[:, step].data
        return np.where(np.logical_and(data >=0, data < self.dz))[0].size

    def _clean_position_matrix(self, step):
        """ Create a new position matrix, discarding trajectories that have left the pore already.
        """

        positions = sparse_matrix(self.positions.shape)
        nonzero = self.positions[:, step].nonzero()[0]

        positions[:nonzero.size, :(positions.shape[1] - step)] = self.positions[nonzero, step:]

        first_free_slot = nonzero.size

        return positions, first_free_slot

    def _count_flux_out(self, step):
        """ Return the number of particles which left the pore at x = L during this step
        """

        return np.where(self.positions[:, step].data >= self.length)[0].size
        #return np.nonzero(self.positions[:, step] >= self.length)[0].size

    def _update_nparticles_in_pore(self, step):
        """ Update the array keeping track of the number of particles in the pore

        :param step: the frame up until which data should analyzed

        :type step: int
        """

        nonzero = self.positions[:, :step].getnnz(axis=0)

        self.nparticles_in_pore = np.concatenate((self.nparticles_in_pore, nonzero))

    def _update_concentration_profile(self, step):

        concentration = np.zeros([step, self.nbins])
        # data = self.positions.T.data

        for t in range(step):
            concentration[t, :] = np.histogram(self.positions[:, t].data, self.nbins, range=(0, self.length), density=False)[0]

        self.concentration = np.concatenate((self.concentration, concentration))

    def plot_number_particles(self, show=False):

        plt.figure()

        plt.plot(self.time_uniform[:self.steps], self.nparticles_in_pore)

        if show:
            plt.show()

    def plot_average_concentration(self, equil=500, show=False, theoretical=True, save=True):
        """ Plot the average concentration profile along the pore

        :param equil: the frame at which the flux can be considered to be equilibrated

        :type equil: int
        """

        plt.figure()
        bar_edges = np.linspace(0, self.length, self.nbins + 1)

        bin_width = bar_edges[1] - bar_edges[0]
        bar_locations = [b + (bin_width / 2) for b in bar_edges[:-1]]
        plt.bar(bar_locations, self.concentration[equil:, :].mean(axis=0), bin_width, align='center')
        plt.xlabel('Distance along pore axis (nm)', fontsize=14)
        plt.ylabel('Concentration', fontsize=14)

        if theoretical:
            # theoretical profile based on measured flux  c = c0 * (1 - x / L)
            x = np.linspace(0, self.length, 1000)
            c = self.pore_concentration * (1 - x / self.length)
            plt.plot(x, c, '--', color='black', lw=2, label='Theoretical Profile')
            plt.legend(fontsize=14)

        plt.tick_params(labelsize=14)
        plt.tight_layout()
        #plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/stochastic_transport/supporting_figures/brownian_conc_profile.pdf')

        # x = np.linspace(0, self.length, 1000)
        # c0 = self.pore_concentration
        # j = self.pore_concentration
        # plt.plot(x, c0 - j * x)

        if show:
            plt.show()

    def animate_concentration_profile(self, bins=25, step=1):

        self.fig, self.ax = plt.subplots()
        nonzero = np.nonzero(self.positions[:, step] >= 0)[0]

        n, self.bins, self.patches = plt.hist(self.positions[nonzero, 0], bins, range=(0, self.length), density=True)

        ani = animation.FuncAnimation(self.fig, self._animate, blit=True,
                                      frames=iter(np.arange(0, self.positions.shape[1], step)), repeat=False)

        # ani.save('/home/bcoscia/brownian_impulse.gif', writer='imagemagick', fps=3)
        plt.show()

    def _animate(self, frame):

        try:
            nonzero = np.nonzero(self.positions[:, frame] >= 0)[0]
            n, _ = np.histogram(self.positions[nonzero, frame], self.bins, range=(0, self.length), density=True)
        except FloatingPointError:  # happens if there is no data
            n = np.zeros(self.bins)

        for rect, h in zip(self.patches, n):
            rect.set_height(h)

        self.ax.set_title('Frame %d' % frame)

        self.ax.relim()
        self.ax.autoscale_view()
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

        return self.patches


if __name__ == "__main__":

    """ NOTES
    - include first step of trajectory (at zero) for flux simulations. But exclude from histogram method
    - script now maintains concentration based on average of all previoius frames. Might want to update to a moving window since might only want to average equilibrium portion with equilibrium concentrations
    - pore concentration fixed to read 0 < x < dz
    - need to make dz independent from number of bins in concentration profile. Can be done by normalizing conc profile
    - larger dz makes it harder to maintain concentration (?)
    """

    L = 10
    ntraj = 500  # this is the number of trajectories that actually make it to the end
    pore_conc = 1
    bins = 25
    steps = 500000
    equil = int(steps/2)  # use 3/4 of the data
    dz = L / bins  # make dz independent of bins. Will just need to change normalization on plot
    print('dz: %.2f' % dz)
    nparticles = 5
    sigma = 1
    dt = 0.001
    print('sigma per step: %.2f' % (sigma * np.sqrt(dt)))
    save = True
    load = False

    nt = 8

    mfpt = Flux(L, ntraj, dt=dt, sigma=sigma, nbins=bins, nt=nt, save=save, load=load)  # higher fluxes require more trajectories
    mfpt.simulate_flux(pore_concentration=pore_conc, dz=dz, steps=steps, measure_flux=True)
    #mfpt.simulate_flux_constant_particles(n=nparticles, dt=dt, steps=steps, measure_flux=True)
    print(mfpt.flux_out[equil:].mean() / dt)

    mfpt.plot_average_concentration(equil=equil)
    mfpt.plot_number_particles(show=True)
    exit()

    plt.hist(mfpt.passage_times, bins=50, range=(0, 10000))

    # for t in range(3):
    #     plt.plot(mfpt.trajectory_hops[t, :, 0], mfpt.trajectory_hops[t, :, 1])
    plt.show()
