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
import numpy as np

# as a function of H, alpha and truncation parameter, t

alphas = np.linspace(1, 2, 10)
H = np.linspace(0.0, 0.5, 11)
trunc = np.linspace(1, 100, 10)
scale = 0.0363
nrealizations = 10

df = pd.DataFrame(index=np.arange(alphas.size*H.size*trunc.size), columns=('H', 'alpha', 't', 'scale', 'max'))

for i, a in enumerate(alphas):
    for j, h in enumerate(H):
        for k, t in enumerate(trunc):

            flm = FLM(h, a, M=4, N=2**12, scale=scale)
            # short traj with lots of realizations faster than single long traj
            flm.generate_realizations(nrealizations, truncate=t)

            print(flm.noise.max())