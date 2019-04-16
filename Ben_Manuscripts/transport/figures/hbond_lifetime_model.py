#!/usr/bin/env python

from LLC_Membranes.timeseries import correlation
import numpy as np
import matplotlib.pyplot as plt

nT = 10000  # number of independent trajectories to generate
nframes = 1028 # number of frames per simulated trajectory

acf = np.zeros(nframes)
nsuccess = np.zeros(nframes)
lengths = np.zeros(nT)
for t in range(nT):
	start = np.random.randint(0, nframes)
	end = np.random.randint(start, nframes)
	traj = np.zeros([nframes])
	traj[start:end][np.random.randint(0, high=(end - start), size=(end - start))] = 1
	plt.plot(traj)
	plt.show()
	exit()
	acf_ = correlation.acf(traj[start:])
	if not np.isnan(acf_[0]):
		acf[:acf_.size] += acf_ 
		nsuccess[:acf_.size] += 1
	lengths[t] = (end - start)

print(np.mean(lengths))
acf /= nsuccess

plt.plot(acf)
plt.show()
