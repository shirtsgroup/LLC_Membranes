#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


alpha = 1 
p = np.array([1, 1, 1])
p = np.array([0.7, 0.3])
data = np.random.dirichlet(alpha*p, size=10000)

fig, ax = plt.subplots(1, p.size, figsize=(12, 5))

for i in range(p.size):
	ax[i].set_title('$B_%d$' %(i + 1), fontsize=16)
	ax[i].hist(data[:, i], bins=50, range=(0, 1))
	ax[i].tick_params(labelsize=14)

plt.show()
