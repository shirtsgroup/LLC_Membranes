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
from LLC_Membranes.timeseries.flm_sim_params import TruncateLevy

# as a function of H, alpha and truncation parameter, t

alphas = np.linspace(1, 2, 10)
H = np.linspace(0.0, 0.5, 11)
trunc = np.linspace(1, 100, 10)
scale = 0.030
nrealizations = 10
load = False
append = True

if not load:

    df = pd.DataFrame(index=np.arange(alphas.size*H.size*trunc.size), columns=('H', 'alpha', 't', 'scale', 'max'))

    for i, a in enumerate(alphas):
        for j, h in enumerate(H):
            for k, t in enumerate(trunc):
                df_ndx = i * (H.size * trunc.size) + j * trunc.size + k
                print('Iteration # %d' % df_ndx)
                flm = FLM(h, a, M=4, N=2**12, scale=scale, truncate=t, correct_truncation=False)  # this is with correction to hurst
                # short traj with lots of realizations faster than single long traj
                flm.generate_realizations(nrealizations, progress=False)

                df.loc[df_ndx] = [h, a, t, scale, flm.noise.max()]
                # df.loc[df_ndx] = [h, a, t, scale, flm.noise.max(axis=1).mean()]

    if append:
        df = df.append(pd.read_pickle('truncate_levy.pl'))  # append to previous pickle

    df.to_pickle('truncate_levy.pl')

else:

    df = pd.read_pickle('truncate_levy.pl')

h_test = 0.35
alpha_test = 1.4
max_value = 1

# trunc = TruncateLevy(data_pickle='truncate_levy.pl')
# t = trunc.interpolate(h_test, alpha_test, max_value, 0.0363)
flm = FLM(h_test, alpha_test, M=4, N=2**12, scale=scale, truncate=max_value)  # this uses TruncateLevy class
flm.generate_realizations(nrealizations)
print(flm.noise.max())

