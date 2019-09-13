#!/usr/bin/env python

""" The fractional_levy_motion module gets the correlation structure wrong. Use this script to generate a database of
hurst parameter inputs / outputs so that a correction can be made to the input hurst paramter to achieve the desired
autocorrelation.
"""

from LLC_Membranes.timeseries.fractional_levy_motion import FLM
from LLC_Membranes.llclib import fitting_functions
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import gamma
from scipy.optimize import curve_fit
import pandas as pd


def hurst_correction(params, a, b, c):
    """ A potential functional form of surface to fit to correction database.

    :param params:
    :param a:
    :param b:
    :param c:

    :return:
    """
    H, alpha = params

    return H + (1 / alpha) - 0.5 + (H - 0.5)*(a + b * alpha) ** c


nrealizations = 100  # Add more realizations to increase the precision of the database grid points
length = 2**12  # length of each realization
load = True

alphas = np.linspace(1, 2, 25)  # range of alphas
H = np.linspace(-0.5, 0.5, 51)  # as input H --> -inf, output H --> 0

if not load:

    df = pd.DataFrame(index=np.arange(alphas.size*H.size), columns=('H', 'alpha', 'h'))

    for j, a in enumerate(alphas):
        for i, h in enumerate(H):

            df_ndx = j * alphas.size + i
            print('Iteration # %d' % df_ndx, end='\r', flush=True)
            df.loc[df_ndx]['H'] = h
            df.loc[df_ndx]['alpha'] = a
            h += ((1 / a) - 0.5)  # recenter Brownian noise

            flm = FLM(h, a, M=4, N=length)
            flm.generate_realizations(nrealizations, progress=False)
            flm.autocorrelation()
            df.loc[df_ndx]['h'] = np.log(2 * flm.acf[:, 1].mean() + 2) / (2 * np.log(2))

    df.to_pickle('hurst_correction.pl')

else:

    df = pd.read_pickle('hurst_correction.pl')

# Test the database
H_data = np.array(df['H'], dtype=float)  # input H
alpha_data = np.array(df['alpha'], dtype=float)  # input alpha
h_data = np.array(df['h'], dtype=float)

H = 0.425
alpha = 1.95833

flm = FLM(H, alpha, M=4, N=length)  # correction happens when this class is initialized
flm.generate_realizations(100)
flm.plot_autocorrelation()
H_est = np.log(2 * flm.acf[:, 1].mean() + 2) / (2 * np.log(2))
print(H_est)
plt.show()
exit()

# get a and b parameters for correction equation
bounds = ([-np.inf, -np.inf, 1], np.inf)
a, b, c = curve_fit(hurst_correction, (H_data, alpha_data), h_data, bounds=bounds)[0]
# plot surface
fig = plt.figure()
ax = fig.gca(projection='3d')

h = np.linspace(0.05, 0.5, 10)
alph = np.linspace(1, 2, 10)
X, Y = np.meshgrid(h, alph)
Z = hurst_correction((X, Y), a, b, c) - (1 / Y - 0.5) #- a
#Z = hurst_correction((X, Y), 2.5, -1.5, 0.5) - (1 / Y - 0.5)

ax.plot_surface(X, Y, Z)
ax.scatter(H_data, alpha_data, h_data)
X, Y = np.meshgrid(H, alphas)
#ax.plot_surface(X, Y, h_data.reshape(alphas.size, h.size))
plt.show()
exit()
# test

