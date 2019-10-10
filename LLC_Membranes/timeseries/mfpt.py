#!/usr/bin/env python

""" Find the mean first passage time (MFPT) of a type of particle
"""

from LLC_Membranes.timeseries.ctrwsim import CTRW
import matplotlib.pyplot as plt
import numpy as np


class MFPT(CTRW):

    def __init__(self, L, n, steps):
        """ Calculate the mean first passage time for a particle by simulating a continuous time random walk of an
        ensemble of particles until they cross some distance threshold.

        :param L: length of 1D path that particle follows (nm).
        :param n: number of particle trajectories to generate
        :param steps: number of steps to take between trajectory evaluations. After the specified number of steps are
        taken, trajectories that have crossed the distance threshold are discontinued. All else are continued.

        :type L: float
        :type n: int
        :type steps: int
        """

        super().__init__(steps, n)

        self.length = L  # how far a particle needs to travel
        self.passage_time = np.zeros(n)
        self.done = np.zeros(n, dtype=bool)
        self.fixed_steps_trajectories()

        # continue running trajectories until all of them have reached +/- L
        self._check_length()

    def _check_length(self):
        """ Find which particles have reached the distance threshold

        :return:
        """

        #TODO: This needs debugging
        print(np.abs(self.z_interpolated).max(axis=1))
        done = np.where(np.abs(self.z_interpolated).max(axis=1) >= self.length)[0]
        self.done[done] = True

        for d in done:
            self.passage_time[d] = np.argmax(self.z_interpolated[d, :] >= self.length)  # converts boolean and finds first True

        print(done)
        print(self.passage_time[done])


if __name__ == "__main__":

    mfpt = MFPT(10, 10, 1000)

    for t in range(10):
        plt.plot(mfpt.z_interpolated[t, :])
    plt.show()
