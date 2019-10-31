#!/usr/bin/env python

""" Find the mean first passage time (MFPT) of a type of particle
"""

from LLC_Membranes.llclib import timeseries
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import tqdm
from multiprocessing import Pool
from scipy.sparse import csr_matrix as sparse_matrix
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings('ignore', "Changing the sparsity structure of a csr_matrix is expensive. lil_matrix is "
                                      "more efficient.")
    warnings.filterwarnings('ignore', "SparseEfficiencyWarning")


class Flux:

    def __init__(self, L, n, dt=1., sigma=1, nbins=25, nt=1, save=True, load=False):
        """ Calculate the mean first passage time for a particle by simulating a continuous time random walk of an
        ensemble of particles until they cross some distance threshold.

        :param L: length of 1D path (a pore for example) that particle follows (nm).
        :param n: number of particle trajectories to generate. It should be at least the number expected to be in a pore
        at steady state, but more is always better.
        :param dt: timestep to use while contsructing trajectories
        :param sigma: width of hop distribution per unit time. The width of the hop distribution for the draws made each
        timestep is sqrt(dt) * sigma
        :param nbins: number of bins in concentration profile
        :param nt: number of threads to use when generating trajectories
        :param save: save trajectories to disk
        :param load: load trajectories saved to disk


        :type L: float
        :type n: int
        :type steps: int
        """

        self.length = L
        self.ntraj = n
        self.sigma = sigma
        self.trajectories = []
        self.time = []
        self.dt = dt

        print('Generating Trajectories...')

        self.generate_trajectories(nt=nt, save=save)

        self.passage_times = None
        self.pore_concentration = 0  # concentration at pore entrance
        #self.pore_cross_sectional_area = np.pi * (pore_radius ** 2)
        self.inlet_volume = 0
        self.dz = 0
        self.steps = 0

        self.nbins = nbins
        self.concentration = np.zeros([0, self.nbins])
        self.nparticles_in_pore = np.array([])

        #self.flux = self.mass_transfer_coefficient * (self.bulk_concentration - self.pore_concentration)
        self.flux_out = None

        self.done = np.zeros(n, dtype=bool)
        self.positions = None
        self.time_uniform = None
        self.pbar = None

        # for animated plotting
        self.patches = None
        self.bins = None
        self.ax = None
        self.fig = None

    def _trajectory_realizations(self, ntraj, sigma):

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

        time = np.ones(len(x))
        crosses_zero = timeseries.switch_points(x > 0)

        trajectories = []
        times = []

        for i, cross in enumerate(crosses_zero[1:]):

            start = crosses_zero[i] + 1
            if i == 0:
                start = 0

            end = cross + 1  # we want to include the last point

            traj = np.abs(x[start:end])
            traj -= traj[0]

            # Recursion!
            if timeseries.switch_points(traj >= 0).size > 2:
                traj_, time_ = self._crossings(traj)
                trajectories += traj_
                times += time_
            else:
                trajectories += [traj]  # + sigma * np.random.normal())
                times += [np.cumsum(time[start:end]) - 1]

        return trajectories, times

    def _trajectory_realizations2(self, ntraj, sigma):

        trajectories = []
        times = []

        n = 0

        while n < ntraj:
            print('\r%d/%d trajectories' % (n, ntraj), end='')
            walk = [0]

            while 0 <= walk[-1] < self.length:

                walk.append(walk[-1] + sigma * np.random.normal())

            walk = np.array(walk)

            trajectories.append(walk)
            times.append(np.arange(walk.size))

            if walk[-1] >= self.length:
                n += 1
                #self.pbar.update(1)

        return trajectories, times

    def generate_trajectories(self, nt=1, save=True, savename='trajectories.npz'):
        """ Generate Brownian trajectories

        :param nt: number of threads
        """

        #self.pbar = tqdm.tqdm(total=self.ntraj)

        sigma = np.sqrt(self.dt) * self.sigma

        pool = Pool(processes=nt)

        trajectories_per_thread = self.ntraj // nt

        arguments = [(trajectories_per_thread, sigma) for _ in range(nt)]

        result = pool.starmap(self._trajectory_realizations2, arguments)

        for thread in range(nt):
            self.trajectories += result[thread][0]
            self.time += result[thread][1]

        if save:
            np.savez_compressed(savename, trajectories=self.trajectories, time=self.time)

    def simulate_flux(self, pore_concentration=1, dt=1, dz=0.1, steps=2000, measure_flux=False):
        """ Simulate a flux into the pore, starting new particles after a specified time
        """

        self.pore_concentration = pore_concentration

        longest = max([len(x) for x in self.trajectories])
        total_trajectories = len(self.trajectories)
        print(total_trajectories)

        # discretize position versus time
        self.positions = sparse_matrix((self.ntraj, longest + steps * self.pore_concentration))

        self.dz = dz
        #self.inlet_volume = self.pore_cross_sectional_area * self.dz
        self.time_uniform = np.arange(self.positions.shape[1]) * dt
        self.steps = steps

        if measure_flux:
            self.flux_out = np.zeros(self.steps)

        traj_no = 0
        start = 0
        naddtot = []

        print('Simulating flux by holding the interface concentration constant at %d particles' % self.pore_concentration)
        for step in tqdm.tqdm(range(steps), unit=' Time Steps'):

            # figure out how many solutes to add in order to keep concentration constant
            nadd = self.pore_concentration - self._get_inlet_concentration(start)

            if nadd < 0:  # pore concentration too high. Don't add any
                nadd = 0

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

            traj_no += nadd
            start += 1

        plt.plot(naddtot)
        #plt.plot(self.flux_out)

        previous_step = self.nparticles_in_pore.size
        self._update_nparticles_in_pore(steps - previous_step)
        self._update_concentration_profile(steps - previous_step)

        #self.animate_concentration_profile()

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
        return np.where(self.positions[:, step].data < self.dz)[0].size

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

        return np.nonzero(self.positions[:, step] >= self.length)[0].size

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
        plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/stochastic_transport/supporting_figures/brownian_conc_profile.pdf')

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

    L = 5
    ntraj = 1000  # this is the number of trajectories that actually make it to the end
    pore_conc = 5
    bins = 25
    steps = 400000
    equil = int(steps/2)  # use 3/4 of the data
    dz = L / bins  # when dz changes, so does the inlet concentration
    nparticles = 100
    sigma = 1
    dt = 0.0005
    save = False
    load = False

    nt = 8

    mfpt = Flux(L, ntraj, dt=dt, sigma=sigma, nbins=bins, nt=nt, save=save, load=load)  # higher fluxes require more trajectories
    mfpt.simulate_flux(pore_concentration=pore_conc, dz=dz, steps=steps)
    #print(mfpt.flux_out[equil:].mean())
    #mfpt.simulate_constant_flux(nparticles=nparticles, dz=dz, steps=steps)
    mfpt.plot_average_concentration(equil=equil)
    mfpt.plot_number_particles(show=True)
    exit()

    plt.hist(mfpt.passage_times, bins=50, range=(0, 10000))

    # for t in range(3):
    #     plt.plot(mfpt.trajectory_hops[t, :, 0], mfpt.trajectory_hops[t, :, 1])
    plt.show()
