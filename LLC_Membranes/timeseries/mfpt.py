#!/usr/bin/env python

""" Find the mean first passage time (MFPT) of a type of particle
"""

from LLC_Membranes.timeseries.ctrwsim import CTRW
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import tqdm


class MFPT(CTRW):

    def __init__(self, L, n, k=0.5, bulk_concentration=4, pore_concentration=4, pore_radius=0.6, nbins=25):
        """ Calculate the mean first passage time for a particle by simulating a continuous time random walk of an
        ensemble of particles until they cross some distance threshold.

        :param L: length of 1D path that particle follows (nm).
        :param n: number of particle trajectories to generate. It should be at least the number expected to be in a pore
        at steady state, but more is always better.
        :param steps: number of steps to take between trajectory evaluations. After the specified number of steps are
        taken, trajectories that have crossed the distance threshold are discontinued. All else are continued.

        :type L: float
        :type n: int
        :type steps: int
        """

        super().__init__(L, n, dwell_dist=None)

        print('Generating Trajectories...')
        self.passage_times = self.mean_first_passage_time()  # Generates trajectories, inherited from CTRW

        self.bulk_concentration = bulk_concentration
        self.mass_transfer_coefficient = k
        self.pore_concentration = pore_concentration  # concentration at pore entrance
        self.pore_cross_sectional_area = np.pi * (pore_radius ** 2)
        self.inlet_volume = 0
        self.dz = 0
        self.steps = 0

        self.nbins = nbins
        self.concentration = np.zeros([0, self.nbins])
        self.nparticles_in_pore = np.array([])

        self.flux = self.mass_transfer_coefficient * (self.bulk_concentration - self.pore_concentration)

        self.done = np.zeros(n, dtype=bool)
        self.positions = None
        self.time_uniform = None

        # for animated plotting
        self.patches = None
        self.bins = None
        self.ax = None
        self.fig = None

    def simulate_flux(self, dt=1, dz=0.1, steps=2000):
        """ Simulate a flux into the pore, starting new particles after a specified time
        """

        # discretize position versus time
        self.positions = np.zeros([self.ntraj, int(np.ceil(np.max(self.passage_times))) +
                                   steps * self.pore_concentration]) - 1  # - 1 to help find inactive trajectories.
        self.dz = dz
        self.inlet_volume = self.pore_cross_sectional_area * self.dz
        self.time_uniform = np.arange(self.positions.shape[1]) * dt
        self.steps = steps

        traj_no = 0
        start = 0
        naddtot = []

        print('Simulating flux by holding the interface concentration constant at %d particles' % self.pore_concentration)
        for step in tqdm.tqdm(range(steps), unit=' Time Steps'):

            # figure out how many solutes to add in order to keep concentration constant

            nadd = self.pore_concentration - self._get_inlet_concentration(start)
            # if step > 0:
            #     print(self._get_inlet_concentration(start - 1))
            #     exit()

            if nadd < 0:  # pore concentration too high. Don't add any
                nadd = 0

            nadd = self.pore_concentration

            naddtot.append(nadd)

            if traj_no + nadd > self.ntraj:

                previous_step = self.nparticles_in_pore.size
                # Record the total number of particles in the pore before modifying position matrix
                self._update_nparticles_in_pore(step - previous_step)

                # Record concentration profile in pore before modifying position matrix
                self._update_concentration_profile(step - previous_step)

                # Remove trajetories that are already finished
                self.positions, traj_no = self._clean_position_matrix(step - previous_step)
                start = 0

                #current_location += previous_step

            for particle in range(traj_no, traj_no + nadd):

                random_trajectory_index = np.random.choice(self.ntraj)  # randomly choose a trajectory

                traj = self.trajectories[random_trajectory_index]
                time = self.time[random_trajectory_index]
                last = np.argmin(np.abs(time[-1] - self.time_uniform))
                tu = self.time_uniform[:(last + 1)]

                # find the indices of the time point closest to the interpolated uniform time series
                time_index = np.argmin(np.abs(time[:, np.newaxis] - tu), axis=0)
                # make sure that jumps don't occur in interpolated trajectory until jumps occur in real trajectory
                time_index[np.where(tu - time[time_index] < 0)[0]] -= 1
                time_index[np.where(time_index < 0)] += 1  # Relax above restriction for first time point

                interpolated = traj[time_index]  # TODO: this has issues for trajectories starting at zero

                self.positions[particle, start:interpolated.size + start] = interpolated
                # print(self._get_inlet_concentration(start))

            traj_no += nadd
            start += 1

        plt.plot(naddtot)

        previous_step = self.nparticles_in_pore.size
        self._update_nparticles_in_pore(steps - previous_step)
        self._update_concentration_profile(steps - previous_step)

        #self.animate_concentration_profile()

    def simulate_constant_flux(self, nparticles=100, dt=1, dz=0.1, steps=2000):
        """ Simulate a flux into the pore, starting new particles after a specified time
        """

        # discretize position versus time
        self.positions = np.zeros([self.ntraj, int(np.ceil(np.max(self.passage_times))) +
                                   steps * self.pore_concentration]) - 1  # - 1 to help find inactive trajectories.
        self.dz = dz
        self.inlet_volume = self.pore_cross_sectional_area * self.dz
        self.time_uniform = np.arange(self.positions.shape[1]) * dt
        self.steps = steps

        traj_no = 0
        start = 0
        naddtot = []
        total_trajectories = len(self.trajectories)
        self.nparticles_in_pore = np.zeros(steps) - 1
        previous_step = 0
        print('Simulating flux by maintaing a constant number of particles in the pore')
        for step in tqdm.tqdm(range(steps), unit=' Time Steps'):

            # figure out how many solutes to add in order to keep concentration constant
            nadd = nparticles - self._get_nparticles_in_pore(start)

            if nadd < 0:  # pore concentration too high. Don't add any
                nadd = 0
            # print(self._get_nparticles_in_pore(start), start, nadd)
            # if step > 0:
            #     positive = np.where(self.positions[:, start] >= 0)[0]
            #     print(positive.size)
            #     exit()

            naddtot.append(nadd)

            if traj_no + nadd > self.ntraj:

                # # Record the total number of particles in the pore before modifying position matrix
                # self._update_nparticles_in_pore(step - previous_step)

                # Record concentration profile in pore before modifying position matrix
                self._update_concentration_profile(step - previous_step)

                # Remove trajetories that are already finished
                self.positions, traj_no = self._clean_position_matrix(step - previous_step)
                start = 0

                previous_step = step

            for particle in range(traj_no, traj_no + nadd):

                random_trajectory_index = np.random.choice(total_trajectories)  # randomly choose a trajectory

                traj = self.trajectories[random_trajectory_index]
                time = self.time[random_trajectory_index]
                last = np.argmin(np.abs(time[-1] - self.time_uniform))
                tu = self.time_uniform[:(last + 1)]

                # find the indices of the time point closest to the interpolated uniform time series
                time_index = np.argmin(np.abs(time[:, np.newaxis] - tu), axis=0)
                # make sure that jumps don't occur in interpolated trajectory until jumps occur in real trajectory
                time_index[np.where(tu - time[time_index] < 0)[0]] -= 1
                time_index[np.where(time_index < 0)] += 1  # Relax above restriction for first time point

                interpolated = traj[time_index]
                # print(time_index[:4])
                # exit()
                interpolated = traj

                self.positions[particle, start:interpolated.size + start] = interpolated
                # print(traj[:2])
                # print(interpolated[:2])
                # print(self.positions[particle, start:start + 3])
                # exit()
            # print(self._get_nparticles_in_pore(start), start)
            # if step > 5:
            #     exit()
            # if start > 1:
            #     negative = np.nonzero(self.positions[:, (start + 1)] < 0)[0]
            #     print(negative.size)
            #     print(self.positions[negative, (start + 1)])
            #     exit()
            #     print(self._get_nparticles_in_pore(start + 1))

            self.nparticles_in_pore[step] = self._get_nparticles_in_pore(start)

            # if start > 1:
            #     print(self.nparticles_in_pore[:(start + 1)])
            #     exit()
            # print(self.positions)

            traj_no += nadd
            start += 1

        plt.plot(naddtot)
        print(previous_step)
        #self._update_nparticles_in_pore(steps - previous_step)
        self._update_concentration_profile(steps - previous_step)

    def _get_inlet_concentration(self, step):
        """ Number of particles in inlet slab

        :param step: step number

        :type step: int

        :return: number of particles in inlet slab
        """

        nonzero = np.nonzero(self.positions[:, step] >= 0)[0]

        return len(np.where(self.positions[nonzero, step] < self.dz)[0])

    def _clean_position_matrix(self, step):
        """ Create a new position matrix, discarding trajectories that have left the pore already.
        """

        positions = np.zeros_like(self.positions) - 1
        nonzero = np.nonzero(self.positions[:, step] >= 0)[0]

        positions[:nonzero.size, :(positions.shape[1] - step)] = self.positions[nonzero, step:]

        first_free_slot = nonzero.size

        return positions, first_free_slot

    def _get_nparticles_in_pore(self, step):

        return np.nonzero(self.positions[:, step] >= 0)[0].size

    def _update_nparticles_in_pore(self, step):
        """ Update the array keeping track of the number of particles in the pore

        :param step: the frame up until which data should analyzed

        :type step: int
        """

        nonzero = np.zeros(step)
        for t in range(step):
            nonzero[t] = self._get_nparticles_in_pore(t)

        self.nparticles_in_pore = np.concatenate((self.nparticles_in_pore, nonzero))

    def _update_concentration_profile(self, step):

        concentration = np.zeros([step, self.nbins])

        for t in range(step):
            nonzero = np.nonzero(self.positions[:, t] >= 0)[0]
            concentration[t, :] = np.histogram(self.positions[nonzero, t], self.nbins, range=(0, self.length),
                                               density=False)[0]
            # concentration[t, :], bins = np.histogram(self.positions[:, t], self.nbins, range=(-self.length, self.length),
            #                                    density=False)

        self.concentration = np.concatenate((self.concentration, concentration))

    def plot_number_particles(self, show=False):

        plt.figure()

        plt.plot(self.time_uniform[:self.steps], self.nparticles_in_pore)

        if show:
            plt.show()

    def plot_average_concentration(self, equil=500, show=False):
        """ Plot the average concentration profile along the pore

        :param equil: the frame at which the flux can be considered to be equilibrated

        :type equil: int
        """

        plt.figure()
        bar_edges = np.linspace(0, self.length, self.nbins + 1)
        bin_width = bar_edges[1] - bar_edges[0]
        bar_locations = [b + bin_width for b in bar_edges[:-1]]
        plt.bar(bar_locations, self.concentration[equil:, :].mean(axis=0), bin_width)
        plt.xlabel('Distance along pore axis (nm)', fontsize=14)
        plt.ylabel('Concentration', fontsize=14)
        plt.tick_params(labelsize=14)

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

    #     # continue running trajectories until all of them have reached +/- L
    #     self._check_length()
    #
    # def _check_length(self):
    #     """ Find which particles have reached the distance threshold
    #     """
    #
    #     done = np.where(np.abs(self.trajectories[..., 1]).max(axis=1) >= self.length)[0]
    #     self.done[done] = True
    #
    #     for d in done:
    #         crossing_time = np.argmax(np.abs(self.trajectories[d, :, 1]) >= self.length)
    #         self.passage_time[d] = self.trajectories[d, crossing_time, 0]  # converts boolean and finds first True
    #
    #     print(done)
    #     print(self.passage_time[done])


if __name__ == "__main__":
    
    L = 10
    ntraj = 5000
    pore_conc = 10
    bins = 200
    steps = 10000
    equil = int(steps/2)
    dz = 0.05
    nparticles = 100

    mfpt = MFPT(L, ntraj, pore_concentration=pore_conc, nbins=bins)  # higher fluxes require more trajectories
    mfpt.simulate_flux(dz=dz, steps=steps)
    #mfpt.simulate_constant_flux(nparticles=nparticles, dz=dz, steps=steps)
    mfpt.plot_average_concentration(equil=equil)
    mfpt.plot_number_particles(show=True)
    exit()

    plt.hist(mfpt.passage_times, bins=50, range=(0, 10000))

    # for t in range(3):
    #     plt.plot(mfpt.trajectory_hops[t, :, 0], mfpt.trajectory_hops[t, :, 1])
    plt.show()
