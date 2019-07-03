#! /usr/bin/env python

from LLC_Membranes.setup import hexagonal_build
from unittest import TestCase as tc
import numpy as np


class TestSetup():

    def test_build(self):

        sys = hexagonal_build.BuildHexagonal('NAcarb11V', 4, 45, 60, 1, tilt=0)  # name, npores, p2p, pore_alpha, pore_radius, tilt=0

        assert sum([sys.tilt, sys.pore_radius]) == 1, 'Should be 1'

        np.testing.assert_almost_equal(168.374582405, np.linalg.norm(sys.pore_centers), decimal=7, err_msg='Should be %.8f' % 168.374582405)

# self.tilt = tilt
# self.xyz = np.zeros([0, 3])
# self.names = []
# self.all_residues = []
# self.pore_radius = pore_radius
#
# # currently only implemented for 4 pores
# self.pore_centers = np.zeros([npores, 2])
#
# pore_alpha_radians = pore_alpha * (np.pi / 180)
#
# self.pore_centers[1, :] = [p2p * np.cos(pore_alpha_radians), p2p * np.sin(pore_alpha_radians)]
# self.pore_centers[2, :] = [p2p * np.cos(pore_alpha_radians) + p2p, p2p * np.sin(pore_alpha_radians)]
# self.pore_centers[3, :] = [p2p, 0]
#
# # center pores in box
# self.pore_centers[:, 0] += (p2p / 2) * (1 + np.cos(pore_alpha_radians))
# self.pore_centers[:, 1] += (p2p / 2) * np.sin(pore_alpha_radians)


if __name__ == "__main__":

    setup = TestSetup()
    setup.test_build()
