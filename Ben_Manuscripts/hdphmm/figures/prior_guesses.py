#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from hdphmm.generate_timeseries import GenARData
from LLC_Membranes.llclib import file_rw
from scipy.stats import norm

com, dt = file_rw.load_object('/home/ben/github/hdphmm/notebooks/comr.pl')

t = np.arange(com.shape[0]) * dt / 1000

traj = com[:, 2, 1]
mu = traj.min() + ((traj.max() - traj.min()) / 2)
sigma = (traj.max() - mu) / 2

x = np.linspace(traj.min() - 0.5, traj.max() + 0.5, 1000)

plt.figure(figsize=(7, 3))

plt.xlabel('Time (ns)', fontsize=14)
plt.ylabel('$z$-coordinate', fontsize=14)
plt.plot(t, traj, lw=2)
plt.plot(1000*norm.pdf(x, mu, sigma), x, color='black', lw=2)
plt.ylim(traj.min() -0.5, traj.max() +0.5)
twosigma = mu + 2*sigma
minustwosigma = mu - 2*sigma
plt.plot([1000*norm.pdf(twosigma, mu, sigma), t[traj.argmax()]], [twosigma, twosigma], '--', color='black') 
plt.plot([1000*norm.pdf(minustwosigma, mu, sigma), t[traj.argmin()]], [minustwosigma, minustwosigma], '--', color='black') 
plt.tight_layout()
plt.savefig('prior_guesses.pdf')
plt.show()

