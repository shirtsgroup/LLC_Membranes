#!/usr/bin/env python

""" When simulating fractional levy motion, truncating the base levy distribution is not trivial. The value you
provide to fractional_levy_motion.py truncates the initial distribution from which random draws are pulled but
correlation structure is added using a fourier transform which changes the maximum drawn values. This script populates
a pandas dataframe with input truncation parameters and the average resulting max value actually observed

"""

import numpy as np
import pandas as pd
from LLC_Membranes.timeseries.fractional_levy_motion import FLM
import matplotlib.pyplot as plt

# as a function of H, alpha and truncation parameter, t

t = 20
flm = FLM(0.4, 1.75, M=4, N=2**12, scale=0.05)
# short traj with lots of realizations faster than single long traj
flm.generate_realizations(10, truncate=t)
flm.plot_marginal(show=True)

print(flm.noise.max())