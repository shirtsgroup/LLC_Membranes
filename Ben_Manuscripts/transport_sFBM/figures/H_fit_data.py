#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# get more data with: for i in $(seq 0 100); do ctrwsim.py -H 0.2 -ntraj 24 -steps 2000 -ensemble -hop fbm -fix_time -acf >> H_estimates.txt; done

dip = []
fit = []

with open('H_estimates.txt', 'r') as f:
	for line in f:
		if line.count('Fit H') == 1:
			fit.append(float(line.split(':')[1]))
		elif line.count('First dip') == 1:
			dip.append(float(line.split(':')[1]))

plt.hist(dip, bins=25, label='Dip')
plt.hist(fit, bins=25, label='Fit')

print('Dip mean : %.3f' % np.mean(dip))
print('Fit mean : %.3f' % np.mean(fit))

plt.show()
