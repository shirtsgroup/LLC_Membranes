#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

centers = [[0, 0], [1, 0]]
r = 0.2

points = np.zeros([10, 2])
thetas = np.linspace(0, 2*np.pi - (2*np.pi/5), 5)

for i in range(5):
	x = r*np.cos(thetas[i])
	y = r*np.sin(thetas[i])
	points[i, :] = [x, y]

thetas += np.random.uniform(0, 2*np.pi)
for i in range(5):
	x = r*np.cos(thetas[i]) + centers[1][0]
	y = r*np.sin(thetas[i]) + centers[1][1]
	points[i + 5, :] = [x, y]

plt.scatter(points[:, 0], points[:, 1])
plt.xlim(-0.5, 1.5)
plt.ylim(-1, 1)
plt.axis('off')
plt.savefig('different_rotation.png')
plt.show()	
